

#weave(joinpath(dirname(pathof(MetidaNCA)), "..", "docs", "validation_report.jmd");
weave(joinpath(dirname(@__FILE__), "validation_report.jmd");
#doctype = "pandoc2pdf",
doctype = "pandoc2pdf",
#out_path = joinpath(dirname(pathof(MetidaNCA)), "..", "docs", "src"),
out_path = joinpath(dirname(@__FILE__), "src"),
pandoc_options=["--toc", "-V colorlinks=true" , "-V linkcolor=blue", "-V urlcolor=red", "-V toccolor=gray"])

rm(joinpath(dirname(@__FILE__), "src", "validation_report.aux"); force=true)
rm(joinpath(dirname(@__FILE__), "src", "validation_report.log"); force=true)
rm(joinpath(dirname(@__FILE__), "src", "validation_report.out"); force=true)

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

=#
