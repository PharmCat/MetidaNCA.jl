################################################################################
# Build the validation report.
################################################################################

include(joinpath(dirname(@__FILE__), "acceptance.jl"))

const VALIDATION_STRICT = lowercase(get(ENV, "VALIDATION_STRICT", "true")) ∉ ("false", "0", "no")

# Start from a clean log, a clean journal and a clean output file, so that
# artefacts of a previous build cannot be mistaken for the result of this one.
acceptance_init()

const VALIDATION_PDF = joinpath(dirname(@__FILE__), "src", "validation_report.pdf")
rm(VALIDATION_PDF; force = true)

weave(joinpath(dirname(@__FILE__), "validation_report.jmd");
doctype = "pandoc2pdf",
out_path = joinpath(dirname(@__FILE__), "src"),
pandoc_options=["--toc", "-V colorlinks=true" , "-V linkcolor=blue", "-V urlcolor=red",
"-V toccolor=gray", "--number-sections"])

# LaTeX intermediates are removed before the gate, so that an aborted build
# leaves the same directory state as a successful one.
for ext in ("aux", "log", "out")
    rm(joinpath(dirname(@__FILE__), "src", "validation_report." * ext); force = true)
end

# The journal is written here as well as by the Conclusion chunk: if that chunk
# errored out, the chunk's own copy never happened, and the journal is the one
# artefact that must exist regardless.
acceptance_write_journal()

# Weave also swallows a failure of the PDF toolchain: it prints the error and
# returns normally, which would leave a build that produced no report at all
# looking successful. This is a toolchain check, not an acceptance check, and it
# is reported separately - hence before the verdict.
if !isfile(VALIDATION_PDF)
    msg = "Validation report was not produced at $VALIDATION_PDF. " *
          "Check the pandoc / LaTeX output above. " *
          "The comparisons that did run are in $ACCEPT_LOG, the divergences in $ACCEPT_JOURNAL."
    VALIDATION_STRICT ? error(msg) : @warn msg
end

# Prints the verdict and the full divergence journal. Raises - and therefore
# fails the documentation build - when any comparison is outside the acceptance
# criteria, when the test suite failed, or when a required scenario did not
# report, unless VALIDATION_STRICT=false.
const VALIDATION_OK = acceptance_verdict(; strict = VALIDATION_STRICT)

#=
using Documenter, MetidaNCA, Weave, PrettyTables, CSV, DataFrames
weave(joinpath(dirname(@__FILE__), "validation_report.jmd");
#doctype = "pandoc2pdf",
doctype = "pandoc2pdf",
out_path = ppath,
pandoc_options=["--toc", "-V colorlinks=true" , "-V linkcolor=blue", "-V urlcolor=red", "-V toccolor=gray"])

mainfont: romanuni.ttf
sansfont: NotoSans-Regular.ttf
monofont: NotoSansMono-Regular.ttf
mathfont: texgyredejavu-math.otf

- name: Install TeXlive and Pandoc
        run: sudo apt-get install pandoc texlive texlive-publishers texlive-science latexmk texlive-xetex texlive-latex-base texlive-fonts-recommended
=#
