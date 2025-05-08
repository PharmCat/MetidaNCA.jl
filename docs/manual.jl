weave(joinpath(dirname(@__FILE__), "manual_ru.jmd");
doctype = "pandoc2pdf",

out_path = joinpath(dirname(@__FILE__), "src"),
pandoc_options=["--toc", "-V colorlinks=true" , "-V linkcolor=blue", "-V urlcolor=red", "-V lang=ru", "-V babel-lang=russian", "-V mainfont=Arial",
"-V toccolor=gray", "--number-sections"])

rm(joinpath(dirname(@__FILE__), "src", "manual_ru.aux"); force=true)
rm(joinpath(dirname(@__FILE__), "src", "manual_ru.log"); force=true)
rm(joinpath(dirname(@__FILE__), "src", "manual_ru.out"); force=true)

#weave(joinpath(dirname(@__FILE__), "manual_ru.jmd"); out_path = "D:\\Temp")