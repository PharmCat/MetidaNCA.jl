

#weave(joinpath(dirname(pathof(MetidaNCA)), "..", "docs", "validation_report.jmd");
weave(joinpath(dirname(@__FILE__), "validation_report.jmd");
#doctype = "pandoc2pdf",
doctype = "github",
#out_path = joinpath(dirname(pathof(MetidaNCA)), "..", "docs", "src"),
out_path = joinpath(dirname(@__FILE__), "src"),
pandoc_options=["--toc", "-V colorlinks=true" , "-V linkcolor=blue", "-V urlcolor=red", "-V toccolor=gray"])
