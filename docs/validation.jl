################################################################################
# Build the validation report.
################################################################################

include(joinpath(dirname(@__FILE__), "acceptance.jl"))

const VALIDATION_STRICT = lowercase(get(ENV, "VALIDATION_STRICT", "true")) ∉ ("false", "0", "no")

const DOCS_DIR  = dirname(@__FILE__)
const SRC_DIR   = joinpath(DOCS_DIR, "src")
# The Markdown intermediate does NOT go into docs/src: Documenter treats every
# .md file there as a page and objects to one that is not listed in `pages`.
const BUILD_DIR = joinpath(DOCS_DIR, "build_md")

const VALIDATION_JMD = joinpath(DOCS_DIR, "validation_report.jmd")
const VALIDATION_MD  = joinpath(BUILD_DIR, "validation_report.md")
const VALIDATION_PDF = joinpath(SRC_DIR, "validation_report.pdf")

const VALIDATION_BUILD_LOG   = joinpath(SRC_DIR, "validation_build.log")
const VALIDATION_PDF_FAILURE = joinpath(SRC_DIR, "validation_pdf_failure.txt")

# ----------------------------------------------------------------- preflight ---
# Packages the report's code chunks need. Weave renders a failing `using` into
# the document and carries on, so one package missing from the active
# environment produces a full-length report built entirely of error messages
# and an acceptance log with nothing in it. Checking first turns that into one
# line, before anything is built. This is an installation-qualification check,
# not an acceptance check, so it is reported separately.
const REPORT_PACKAGES = ["Dates", "DataFrames", "CSV", "PrettyTables", "Latexify",
                         "Statistics", "Plots", "MetidaNCA", "Pkg"]

function validation_preflight(; strict = true)
    absent = filter(p -> isnothing(Base.find_package(p)), REPORT_PACKAGES)
    isempty(absent) && return true
    msg = """
    VALIDATION CANNOT RUN: the report needs $(join(absent, ", ")), not available in the
    active project ($(Base.active_project())).

    Every code chunk of the report would fail, the acceptance log would stay empty, and the
    resulting PDF would consist of error messages.

    Instantiate the documentation environment and build with it:

        julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
        julia --project=docs docs/make.jl
    """
    strict ? error(msg) : @warn msg
    return false
end

"""
    validation_subject_check(; strict = true) -> Bool

Check that the MetidaNCA the report will exercise is the one in this repository.

This is the first question a validation report has to answer - *which* software
was qualified - and it is easy to get wrong without noticing. `docs/Project.toml`
declares MetidaNCA as an ordinary dependency, so unless the documentation
environment has `Pkg.develop(path = ".")` applied to it, the environment resolves
MetidaNCA to a copy from the registry. The build then succeeds and produces a
report that qualifies a released version while appearing to qualify the working
tree, and every parameter in it belongs to code the repository does not contain.

The giveaway in the log is a stack frame under `.julia/packages/MetidaNCA/`
rather than under the repository.
"""
function validation_subject_check(; strict = true)
    pkg = Base.find_package("MetidaNCA")
    isnothing(pkg) && return false          # already reported by the preflight
    pkgroot  = normpath(abspath(joinpath(dirname(pkg), "..")))
    reporoot = normpath(abspath(joinpath(DOCS_DIR, "..")))
    # Same directory, up to the separator and case conventions of the platform.
    same = lowercase(rstrip(pkgroot, ['/', '\\'])) == lowercase(rstrip(reporoot, ['/', '\\']))
    same && return true
    msg = """
    VALIDATION WOULD QUALIFY THE WRONG CODE.

    The report is in   $reporoot
    but MetidaNCA resolves to
                       $pkgroot

    The documentation environment is using a copy of the package from the registry, not this
    repository, so the report would describe a released version while appearing to describe
    the working tree. It would also read its data and reference files out of that copy, which
    may not contain files added since the release.

    Point the documentation environment at this repository:

        julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

    and build again. To qualify a released version deliberately, check that version out and
    build from its own working tree.
    """
    strict ? error(msg) : @warn msg
    return false
end

validation_preflight(; strict = VALIDATION_STRICT)
validation_subject_check(; strict = VALIDATION_STRICT)

# ------------------------------------------------------------- step 1: Weave ---

# Start from a clean log, a clean journal and clean output, so that artefacts of
# a previous build cannot be mistaken for the result of this one.
acceptance_init()
rm(VALIDATION_PDF; force = true)
rm(VALIDATION_PDF_FAILURE; force = true)
rm(BUILD_DIR; force = true, recursive = true)
mkpath(BUILD_DIR)

# Weave's own messages, and anything a chunk writes to stderr, are kept on disk:
# a warning raised while a table was being built is otherwise lost in the
# scrollback of a build that takes minutes.
open(VALIDATION_BUILD_LOG, "w") do logio
    redirect_stderr(logio) do
        weave(VALIDATION_JMD; doctype = "pandoc", out_path = BUILD_DIR)
    end
end

"""
    build_log_tail(n = 40) -> String

The last `n` lines of the Weave log that carry a message, with the chunk
progress chatter dropped.
"""
function build_log_tail(n::Int = 40)
    isfile(VALIDATION_BUILD_LOG) || return "(no build log at $VALIDATION_BUILD_LOG)"
    lines = filter(readlines(VALIDATION_BUILD_LOG)) do l
        s = strip(l)
        !isempty(s) && !occursin("Weaving chunk", s) && !occursin("progress =", s) &&
            !occursin("Weaved all chunks", s)
    end
    isempty(lines) && return "(build log holds nothing but progress messages)"
    return join(lines[max(1, end - n + 1):end], "\n")
end

# The journal is written here as well as by the Conclusion chunk: if that chunk
# errored out, the chunk's own copy never happened, and the journal is the one
# artefact that must exist regardless.
acceptance_write_journal()

# ------------------------------------------------------------ step 2: pandoc ---
# The same command Weave's Pandoc2PDF format issues, with two differences: the
# output is captured instead of discarded, and the `-V` options are passed as
# separate arguments. Weave is given them as single strings with a space inside
# ("-V colorlinks=true"), which reaches pandoc as one argument and sets a
# variable whose name begins with a space - that is, it never worked.

# Weave's own header template is a single \usepackage{geometry} line; it is
# reproduced here so the conversion does not reach into Weave's internals and
# can be repeated by hand.
const PANDOC_HEADER = joinpath(DOCS_DIR, "pandoc_header.tex")
isfile(PANDOC_HEADER) || write(PANDOC_HEADER,
    "\\usepackage[a4paper,text={16.5cm,25.2cm},centering]{geometry}\n")

const PANDOC_OPTIONS = ["--toc", "--number-sections",
                        "-V", "colorlinks=true",
                        "-V", "linkcolor=blue",
                        "-V", "urlcolor=red",
                        "-V", "toccolor=gray"]

"""
    validation_pandoc(md, pdf) -> (ok::Bool, output::String)

Convert the rendered Markdown to PDF and return whether a PDF came out,
together with everything pandoc and the TeX engine had to say. pandoc runs with
the Markdown's directory as the working directory, because the figures are
referenced relative to it.
"""
function validation_pandoc(md, pdf; options = PANDOC_OPTIONS, header = PANDOC_HEADER)
    out = IOBuffer()
    cd(dirname(abspath(md))) do
        cmd = `pandoc -f markdown+raw_tex -s --pdf-engine=xelatex --highlight-style=tango
               $options --include-in-header=$(abspath(header))
               -V fontsize=12pt -o $(abspath(pdf)) $(basename(md))`
        println(out, "command: ", cmd)
        println(out, "working directory: ", pwd())
        println(out)
        try
            run(pipeline(ignorestatus(cmd); stdout = out, stderr = out))
        catch err
            println(out, "\npandoc could not be run at all: ", sprint(showerror, err))
        end
    end
    return isfile(pdf), String(take!(out))
end

const PDF_OK, PANDOC_OUTPUT = validation_pandoc(VALIDATION_MD, VALIDATION_PDF)

# LaTeX intermediates, wherever the engine chose to leave them.
for dir in (BUILD_DIR, SRC_DIR), ext in ("aux", "log", "out", "toc")
    rm(joinpath(dir, "validation_report." * ext); force = true)
end

if !PDF_OK
    figdir = joinpath(BUILD_DIR, "figures")
    figs = isdir(figdir) ? readdir(figdir) : String[]
    report = """
    The validation report was not converted to PDF.

    The code chunks ran. Their comparisons are in
        $ACCEPT_LOG
    and the divergences in
        $ACCEPT_JOURNAL
    What failed is the conversion, after all the computation. The rendered Markdown is kept at
        $VALIDATION_MD
    and the figures at
        $figdir
    ($(length(figs)) file(s)$(isempty(figs) ? " - none, which is itself the problem if the report has figures" : ""))

    ================================ pandoc output =================================
    $PANDOC_OUTPUT
    ================================ end of output =================================

    Weave log ($VALIDATION_BUILD_LOG), last messages:
    --------------------------------------------------------------------------------
    $(build_log_tail(25))
    --------------------------------------------------------------------------------

    Environment:
      julia   $(VERSION) on $(Sys.MACHINE)
      pandoc  $(try first(eachline(IOBuffer(read(`pandoc --version`, String)))) catch; "not callable" end)
      xelatex $(try first(eachline(IOBuffer(read(`xelatex --version`, String)))) catch; "not callable" end)

    To retry the conversion alone, without recomputing anything:

        cd $(BUILD_DIR)
        pandoc -f markdown+raw_tex -s --pdf-engine=xelatex --highlight-style=tango \\
          --toc --number-sections --include-in-header=$PANDOC_HEADER \\
          -V fontsize=12pt -o report.pdf validation_report.md
    """
    write(VALIDATION_PDF_FAILURE, report)
    msg = report * "\n\nThis text is also in $VALIDATION_PDF_FAILURE."
    VALIDATION_STRICT ? error(msg) : @warn msg
end

# ------------------------------------------------------------------- the gate ---
# Prints the verdict and the full divergence journal. Raises - and therefore
# fails the documentation build - when any comparison is outside the acceptance
# criteria, when the test suite failed, or when a required scenario did not
# report, unless VALIDATION_STRICT=false.
const VALIDATION_OK = acceptance_verdict(; strict = VALIDATION_STRICT)

#=
mainfont: romanuni.ttf
sansfont: NotoSans-Regular.ttf
monofont: NotoSansMono-Regular.ttf
mathfont: texgyredejavu-math.otf

- name: Install TeXlive and Pandoc
        run: sudo apt-get install pandoc texlive texlive-publishers texlive-science latexmk texlive-xetex texlive-latex-base texlive-fonts-recommended
=#
