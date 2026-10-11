# =============================================================================
#  MetidaNCA.jl 0.7.1 — тесты, воспроизводящие проблемы из отчёта
# =============================================================================
#  Каждый @testset соответствует номеру проблемы в отчёте.
#  Утверждения записаны так, как ДОЛЖНО быть → на текущем коде они падают.
#  Вспомогательные функции печатают «ожидается / фактически».
#
#  Запуск:
#      julia --project=. metidanca_bugtests.jl
#  или из REPL:
#      include("metidanca_bugtests.jl")
#
#  Требуется: MetidaNCA, CSV, DataFrames, Plots, Unitful, Test
# =============================================================================
import Pkg
path = dirname(@__FILE__)
cd(path)


using Test, CSV, DataFrames, Plots, Unitful
import Main.MetidaNCA: pkimport, upkimport, pdimport, nca, nca!, DoseTime, ElimRange,
                 LimitRule, NoPageSort, applylimitrule!, setkelrange!, setkelauto!,
                 setdosetime!, getkelrange, getkelauto, getkeldata, getdosetime,
                 timefilter, auc_sparse, pkplot, NCAResult

# ---- где лежат csv из репозитория ------------------------------------------
const CSVDIR = let d = get(ENV, "METIDANCA_CSV", "")
    isempty(d) ? normpath(joinpath(path, "csv")) : d
end
@assert isdir(CSVDIR) "не найден каталог с csv: $CSVDIR (задайте ENV[\"METIDANCA_CSV\"])"

pkdata2   = CSV.File(joinpath(CSVDIR, "pkdata2.csv"))   |> DataFrame
missingpk = CSV.File(joinpath(CSVDIR, "missingpk.csv")) |> DataFrame
upkdata   = CSV.File(joinpath(CSVDIR, "upkdata.csv"))   |> DataFrame
pddata    = CSV.File(joinpath(CSVDIR, "pddata.csv"))    |> DataFrame

pk1 = pkdata2[pkdata2.Subject .== 1, :]          # 16 точек, 0 … 72 ч
devnull_io = IOBuffer()

# ---- утилиты вывода ---------------------------------------------------------
const SUMMARY = Tuple{String,Bool}[]

_short(e) = first(split(sprint(showerror, e), '\n'))

"""Сравнить ожидаемое и фактическое, напечатать строку, вернуть Bool."""
function cmpshow(label, expected, actual; atol = 1e-6)
    ok = try
        isa(expected, Number) && isa(actual, Number) ?
            isapprox(expected, actual; atol = atol) : expected == actual
    catch
        false
    end
    printstyled(rpad(label, 52); color = :normal)
    printstyled(ok ? "OK   " : "FAIL "; color = ok ? :green : :red)
    println("ожидается: ", expected, "  |  фактически: ", actual)
    push!(SUMMARY, (label, ok))
    ok
end

"""Выполнить f(); ожидается отсутствие исключения. Печатает исключение, если оно есть."""
function nothrows(label, f)
    ok = true
    try
        f()
        printstyled(rpad(label, 52); color = :normal); printstyled("OK   \n"; color = :green)
    catch e
        ok = false
        printstyled(rpad(label, 52); color = :normal); printstyled("FAIL "; color = :red)
        println("исключение: ", _short(e))
    end
    push!(SUMMARY, (label, ok))
    ok
end

"""Значение по ключу результата или :absent, если ключа нет."""
key(res, k) = haskey(res.result, k) ? res.result[k] : :absent

# Тесты, способные уронить сессию через StackOverflowError, по умолчанию выключены
const RUN_CRASH_TESTS = get(ENV, "METIDANCA_CRASH_TESTS", "0") == "1"

println("\n", "="^100)
println("MetidaNCA — воспроизведение проблем отчёта")
println("="^100, "\n")


# =============================================================================
# №1  Ось времени сдвигается только при tau > 0; точка Cdose вставляется в t = 0
# =============================================================================
@testset "№1  dosetime > 0 при tau = 0 игнорируется в AUC/AUMC" begin
    println("\n--- №1 -------------------------------------------------------------")

    r0    = nca(pk1, :Time, :Concentration; dosetime = DoseTime(dose = 100, time = 0.0), verbose = 3)
    r025  = nca(pk1, :Time, :Concentration; dosetime = DoseTime(dose = 100, time = 0.25), verbose = 3)
    r2    = nca(pk1, :Time, :Concentration; dosetime = DoseTime(dose = 100, time = 2.0), verbose = 3)

    # (a) сдвиг времени дозы обязан уменьшить AUClast: интервал [0; 0.25] не должен учитываться.
    #     Эталон, посчитанный трапециями от t = 0.25 (C(0.25) = 89.4745): 9574.2374375
    @test cmpshow("№1a AUClast при dosetime = 0.25", 9574.2374375, r025[:AUClast]; atol = 1e-4)

    # (b) прямое следствие: результат не должен совпадать с расчётом от t = 0
    @test cmpshow("№1b AUClast(0.25) != AUClast(0)", true, !isapprox(r0[:AUClast], r025[:AUClast]))

    # (c) время дозы 2.0 ч. Эталон трапециями по точкам 2 … 72: 9283.05875
    #     Фактически добавляется фантомная трапеция [0; 2] с C(0) = 0 → +139.962
    @test cmpshow("№1c AUClast при dosetime = 2.0", 9283.05875, r2[:AUClast]; atol = 1e-4)

    # (d) AUMC обязан считаться от времени дозы (docs/parameters.md: «Dose time is
    #     the starting point»). Эталон по (t - 2): 314682.51525
    @test cmpshow("№1d AUMClast при dosetime = 2.0", 314682.51525, r2[:AUMClast]; atol = 1e-2)

    # (e) частичная AUC при ненулевом времени дозы: интервал сдвигается на dosetime,
    #     а сама ось времени — нет, поэтому считается не тот участок.
    rp = nca(pk1, :Time, :Concentration; dosetime = DoseTime(dose = 100, time = 2.0),
             partials = [(2.0, 4.0)], prtext = :last)
    #     ожидается AUC на [2; 4] = 67.388 + 65.23975 + 141.1115 = 273.73925
    #     фактически интервал сдвигается в [0; 2], где стоит фантомная точка → 139.962
    @test cmpshow("№1e partial AUC [2;4] при dosetime = 2.0", 273.73925,
                  key(rp, Symbol("AUC_2.0_4.0")); atol = 1e-3)
end


# =============================================================================
# №2  Partial AUC = NaN, если конец интервала точно совпал с точкой наблюдения
# =============================================================================
@testset "№2  partial AUC → NaN при попадании в точку отбора" begin
    println("\n--- №2 -------------------------------------------------------------")

    subj = pkimport(pk1, :Time, :Concentration; dosetime = DoseTime(dose = 100, time = 0.0))
    r    = nca!(subj; partials = [(1, 3), (1, 7)])      # 3.0 — точка отбора, 7.0 — нет

    # AUC[1;3] линейными трапециями = 88.949 + 76.22225 + 67.388 + 65.23975 = 297.799
    # фактически: конец интервала (3.0) совпал с точкой отбора → ветка else → NaN
    @test cmpshow("№2a AUC_1_3 конечен", true, isfinite(key(r, :AUC_1_3)))
    @test cmpshow("№2b AUC_1_3 значение", 297.799, key(r, :AUC_1_3); atol = 1e-3)
    # контроль: интервал, конец которого лежит МЕЖДУ точками, считается нормально
    @test cmpshow("№2c AUC_1_7 конечен (контроль)", true, isfinite(key(r, :AUC_1_7)))
end


# =============================================================================
# №3  nca! мутирует входной субъект; индексы попадают в чужую систему отсчёта
# =============================================================================
@testset "№3  nca! изменяет kelrange входного субъекта" begin
    println("\n--- №3 -------------------------------------------------------------")

    # (a) сам факт мутации: анализ не должен менять объект-вход
    s1 = pkimport(missingpk, :Time, :Concentration)
    @assert isempty(getkelrange(s1).kelexcl)
    nca!(s1)
    @test cmpshow("№3a kelexcl входного субъекта не изменён", Int[], getkelrange(s1).kelexcl)

    # (b) система отсчёта: einds — индексы ОТФИЛЬТРОВАННОГО вектора,
    #     kelexcl — индексы ИСХОДНОГО. При dosetime > 0 они расходятся.
    #     Интерполировались точки t = 11.0 (missing) и t = 13.0 (NaN).
    s2 = pkimport(missingpk, :Time, :Concentration;
                  dosetime = DoseTime(dose = 100, time = 0.5))
    nca!(s2)
    excl = getkelrange(s2).kelexcl
    got  = isempty(excl) ? Float64[] : Float64.(MetidaNCA.gettime(s2)[excl])
    @test cmpshow("№3b времена исключённых точек", [11.0, 13.0], got)

    # (c) ElimRange{:time}: kelexcl хранит ВРЕМЕНА, а дописываются целые индексы
    s3 = pkimport(missingpk, :Time, :Concentration;
                  elimrange = ElimRange(10.0, 72.0, Float64[]; time = true))
    nca!(s3)
    @test cmpshow("№3c kelexcl (time-режим) не засорён индексами",
                  Float64[], getkelrange(s3).kelexcl)

    # (d) неидемпотентность: kelauto пользователя молча переписывается
    s4 = pkimport(missingpk, :Time, :Concentration;
                  kelauto = false, elimrange = ElimRange(kelstart = 13, kelend = 18))
    nca!(s4)
    @test cmpshow("№3d kelauto входного субъекта сохранён", false, getkelauto(s4))
end


# =============================================================================
# №4  timefilter фильтрует по элементам диапазона, а не по интервалу
# =============================================================================
@testset "№4  timefilter отбрасывает почти все точки" begin
    println("\n--- №4 -------------------------------------------------------------")

    subj = pkimport(pk1, :Time, :Concentration)
    # времена pk1: 0, .5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 24, 48, 72
    f1 = timefilter(subj, (0.5, 12.0))
    # в интервале [0.5; 12] — 12 точек
    @test cmpshow("№4a timefilter((0.5, 12.0)) — число точек", 12, length(f1.time))
    @test cmpshow("№4b оставшиеся времена", [0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0,10.0,12.0],
                  collect(f1.time))

    # тот же эффект для AbstractRange: 0.5/1.5/2.5 не являются элементами 0:1:12
    f2 = timefilter(subj, 0:1:12)
    @test cmpshow("№4c timefilter(0:1:12) — число точек", 12, length(f2.time))
end


# =============================================================================
# №5  step_3_elim!: бесконечная рекурсия при ненайденном kelstart/kelend
# =============================================================================
@testset "№5  ненайденный kelstart → StackOverflowError" begin
    println("\n--- №5 -------------------------------------------------------------")

    t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    c = Union{Missing,Float64}[0.0, 10.0, 8.0, 5.0, missing, missing]
    # хвост из missing обрезается → time_cp = [0,1,2,3]; kelstart = 5 указывает на t = 4.0,
    # которого в отфильтрованном векторе нет → getstimep вернёт nothing
    if RUN_CRASH_TESTS
        subj = pkimport(t, c; kelauto = false, elimrange = ElimRange(kelstart = 5, kelend = 6))
        @test nothrows("№5  nca! при kelstart вне окна", () -> nca!(subj))
    else
        println("№5  пропущен (может уронить сессию). Включить: ENV[\"METIDANCA_CRASH_TESTS\"]=\"1\"")
        println("     Ожидается: понятная ошибка или автоматический переход к kelauto.")
        println("     Фактически: @warn + рекурсивный вызов с тем же kelauto=false → StackOverflowError,")
        println("                 плюс результат рекурсии отбрасывается и вторая проверка тестирует stimep.")
    end
end


# =============================================================================
# №6  Мультиобсервационный субъект рассинхронизируется правилами/фильтрами
# =============================================================================
@testset "№6  LimitRule/drop* ломают мультиобсервационный субъект" begin
    println("\n--- №6 -------------------------------------------------------------")

    df = deepcopy(pk1)
    df.Conc2 = copy(df.Concentration)
    subj = pkimport(df, :Time, [:Concentration, :Conc2])

    # LimitRule(rm = true) удаляет точки из time и из ПЕРВОГО вектора наблюдений
    applylimitrule!(subj, LimitRule(lloq = 130.0, btmax = NaN, atmax = NaN, rm = true))

    @test cmpshow("№6a длина time == длине Conc2",
                  length(subj.time), length(subj.obs.Conc2))
    @test cmpshow("№6b длина Concentration == длине Conc2",
                  length(subj.obs.Concentration), length(subj.obs.Conc2))
    # как следствие — расчёт по второму наблюдению падает на проверке длин
    @test nothrows("№6c nca! по второму наблюдению", () -> nca!(subj, :Conc2))

    # тот же дефект в dropfunction! (dropnanormissing!/dropnan!/dropmissing!)
    df2 = deepcopy(missingpk)
    df2.Conc2 = copy(df2.Concentration)
    subj2 = pkimport(df2, :Time, [:Concentration, :Conc2])
    MetidaNCA.dropnanormissing!(subj2)
    @test cmpshow("№6d dropnanormissing! сохраняет согласованность",
                  length(subj2.time), length(subj2.obs.Conc2))
end


# =============================================================================
# №7  pkimport без sort не сортирует по времени
# =============================================================================
@testset "№7  pkimport(data, time, obs) не сортирует наблюдения" begin
    println("\n--- №7 -------------------------------------------------------------")

    pk1rev = pk1[end:-1:1, :]

    a = nca(pk1,    :Time, :Concentration)[:AUClast]                 # уже отсортировано
    b = nca(pk1rev, :Time, :Concentration)[:AUClast]                 # метод без sort
    c = nca(pk1rev, :Time, :Concentration, :Subject)[1][:AUClast]    # метод с sort — сортирует

    @test cmpshow("№7a AUClast не зависит от порядка строк", a, b; atol = 1e-4)
    @test cmpshow("№7b контроль: метод с sort устойчив", a, c; atol = 1e-4)

    # побочно: pkimport(time, obs) не копирует вектор наблюдений (алиасинг)
    cc     = [0.0, 10.0, 8.0, 5.0]
    cc_ref = copy(cc)
    sj     = pkimport([0.0, 1.0, 2.0, 3.0], cc)
    applylimitrule!(sj, LimitRule(lloq = 6.0, btmax = 0.0, atmax = 0.0))
    @test cmpshow("№7c исходный вектор наблюдений не изменён", cc_ref, cc)
end


# =============================================================================
# №8  Unitful: жёсткие Float64 и литералы
# =============================================================================
@testset "№8  Unitful-ветки" begin
    println("\n--- №8 -------------------------------------------------------------")

    upk = deepcopy(pk1)
    upk.Time          = upk.Time .* u"hr"
    upk.Concentration = upk.Concentration .* u"ng/ml"

    # (a) kelauto = false → findall(x -> x <= 0, obs_cp) сравнивает Quantity с 0
    s_a = pkimport(upk, :Time, :Concentration;
                   kelauto = false, elimrange = ElimRange(kelstart = 10, kelend = 16))
    @test nothrows("№8a nca! с kelauto = false", () -> nca!(s_a))

    # (b) verbose > 0 → pushfirst!(aucpartl, 0.0) в вектор Quantity
    s_b = pkimport(upk, :Time, :Concentration)
    @test nothrows("№8b nca! с verbose = 1", () -> nca!(s_b; verbose = 1, io = IOBuffer()))

    # (c) экстраполяция Ctau: exp(...) без * oneunit(O) → потеря единиц
    s_c = pkimport(upk, :Time, :Concentration;
                   dosetime = DoseTime(dose = 100u"mg", time = 0.0u"hr", tau = 100.0u"hr"))
    @test nothrows("№8c nca! с tau > последнего времени", () -> nca!(s_c))

    # (d) экстраполяция partials: тот же exp(...) без единиц
    s_d = pkimport(upk, :Time, :Concentration;
                   dosetime = DoseTime(dose = 100u"mg", time = 0.0u"hr"))
    @test nothrows("№8d partials с prtext = :extr",
                   () -> nca!(s_d; partials = [(0.0u"hr", 100.0u"hr")], prtext = :extr))

    # (e) auc_sparse: wts = Vector{Float64}
    @test nothrows("№8e auc_sparse с единицами",
                   () -> auc_sparse([0.0, 1.0, 2.0] .* u"hr", [1.0, 2.0, 1.0] .* u"ng/ml"))

    # (f) LimitRule хранит Float64 → сравнение с размерной концентрацией
    s_f = pkimport(upk, :Time, :Concentration)
    @test nothrows("№8f applylimitrule! на Unitful-субъекте",
                   () -> applylimitrule!(s_f, LimitRule(lloq = 1.0, btmax = 0.0, atmax = 0.0)))
end


# =============================================================================
# №9  vpcdata возвращает границы интервала перепутанными
# =============================================================================
@testset "№9  vpcdata: lb и ub поменяны местами" begin
    println("\n--- №9 -------------------------------------------------------------")

    ds = pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    k, means, lb, ub = MetidaNCA.vpcdata(ds)

    @test cmpshow("№9a lb <= means всюду", true, all(lb .<= means))
    @test cmpshow("№9b means <= ub всюду", true, all(means .<= ub))
    @test cmpshow("№9c lb <= ub всюду",    true, all(lb .<= ub))
end


# =============================================================================
# №10–12, 46  Графики
# =============================================================================
@testset "№10-12,46  pkplot" begin
    println("\n--- №10-12, 46 -----------------------------------------------------")

    ds   = pkimport(pkdata2, :Time, :Concentration, [:Subject, :Formulation])
    subj = ds[1]

    # №10 drawdt для PKSubject: _pddtplot! обращается к .time у Vector{DoseTime}
    @test nothrows("№10 pkplot(PKSubject; drawdt = true)",
                   () -> pkplot(subj; drawdt = true))

    # №11 uylims: в pkplot используется необъявленная переменная obsname
    @test nothrows("№11 pkplot(DataSet; uylims = true)",
                   () -> pkplot(ds; uylims = true))

    # №12 _subjplot: yscale вне списка логарифмов → time/obs не присваиваются
    @test nothrows("№12 pkplot(subj; yscale = :identity)",
                   () -> pkplot(subj; yscale = :identity))

    # №46 pkplot(DataSet{NCAResult}) предполагает DataSet{PKPlot}
    dsnca = nca!(ds)
    @test nothrows("№46 pkplot(результаты; onlyplots = true, elim = true)",
                   () -> pkplot(dsnca; onlyplots = true, elim = true))
end


# =============================================================================
# №13–14  Фасад nca(...) теряет/не пропускает keyword-аргументы
# =============================================================================
@testset "№13-14  nca(...) для :ur и :pd" begin
    println("\n--- №13-14 ---------------------------------------------------------")

    # №13 у upkimport/pdimport нет kwargs... — любой лишний keyword ломает вызов
    @test nothrows("№13a nca(:ur) с calcm",
                   () -> nca(upkdata, :st, :et, :conc, :vol; type = :ur, calcm = :luld))
    @test nothrows("№13b nca(:pd) с verbose",
                   () -> nca(pddata, :time, :obs, :subj; type = :pd, bl = 1.5, th = 5.0,
                             verbose = 1, io = IOBuffer()))

    # №14 dosetime — явный keyword nca(), но в pdimport он не передаётся
    r_pd = nca(pddata, :time, :obs, :subj; type = :pd, bl = 1.5, th = 5.0,
               dosetime = DoseTime(dose = 1, time = 0.0, tau = 5.0))
    @test cmpshow("№14a τ-параметры ФД посчитаны", true, haskey(r_pd[1].result, :AUCABLtau))

    # контроль: через pdimport напрямую они появляются
    pd_ok = pdimport(pddata, :time, :obs, :subj; bl = 1.5, th = 5.0,
                     dosetime = DoseTime(dose = 1, time = 0.0, tau = 5.0))
    r_ok  = nca!(pd_ok)
    @test cmpshow("№14b контроль: pdimport(dosetime=...) работает",
                  true, haskey(r_ok[1].result, :AUCABLtau))
end


# =============================================================================
# №15–18  Urine PK
# =============================================================================
@testset "№15-18  UPKSubject" begin
    println("\n--- №15-18 ---------------------------------------------------------")

    # №15 AR считается по НЕотфильтрованным данным
    #     upkdata: conc*vol = 1, 4, 6, 3, 2. При времени дозы 2 ч учитываются
    #     только интервалы (2,6), (6,12), (12,18) → AR = 6+3+2 = 11, Vol = 8
    u15 = upkimport(upkdata, :st, :et, :conc, :vol, :subj;
                    dosetime = DoseTime(dose = 100, time = 2.0))
    r15 = nca!(u15[1])
    @test cmpshow("№15a AR по данным после времени дозы", 11.0, r15[:AR])
    @test cmpshow("№15b Vol (контроль, фильтруется)",      8.0,  r15[:Vol])

    # №16 kelexcl/kelrange для мочи: gettime — вектор кортежей
    @test nothrows("№16a nca! мочи с kelauto = false",
                   () -> nca!(upkimport(upkdata, :st, :et, :conc, :vol, :subj;
                                        kelauto = false,
                                        elimrange = ElimRange(kelstart = 2, kelend = 5))[1]))
    r16a = nca!(upkimport(upkdata, :st, :et, :conc, :vol, :subj;
                          dosetime = DoseTime(dose = 100, time = 0.0))[1])
    r16b = nca!(upkimport(upkdata, :st, :et, :conc, :vol, :subj;
                          dosetime = DoseTime(dose = 100, time = 0.0),
                          elimrange = ElimRange(kelexcl = Int[4]))[1])
    @test cmpshow("№16b исключение точки меняет результат", true,
                  key(r16a, :AUCall) != key(r16b, :AUCall) || key(r16a, :Kel) != key(r16b, :Kel))

    # №17 нулевая точка ставится в t = 0, а не во время дозы
    #     dosetime = 1.5: mtime должно начинаться с 1.5 → AUCall = 1.875 + 5 + 2.5 = 9.375
    #     фактически mtime[1] = 0 → AUCall = 3 + 5 + 2.5 = 10.5
    u17 = upkimport(upkdata, :st, :et, :conc, :vol, :subj;
                    dosetime = DoseTime(dose = 100, time = 1.5))
    r17 = nca!(u17[1])
    @test cmpshow("№17 AUCall при dosetime = 1.5", 9.375, r17[:AUCall]; atol = 1e-6)

    # №18 Rlast/AUClast при нулевой последней скорости выведения
    u18 = upkimport([0.0, 1.0, 2.0], [1.0, 2.0, 3.0], [1.0, 2.0, 0.0], [1.0, 1.0, 1.0])
    r18 = nca!(u18)
    @test cmpshow("№18a Rlast — последняя измеримая скорость", 2.0, r18[:Rlast])
    @test cmpshow("№18b AUClast < AUCall", true, r18[:AUClast] < r18[:AUCall])
end


# =============================================================================
# №19  LimitRule: BLQ определяется как <=, а не <
# =============================================================================
@testset "№19  значение, равное LLOQ, зануляется" begin
    println("\n--- №19 ------------------------------------------------------------")

    subj = pkimport([0.0, 1.0, 2.0, 3.0], [0.0, 10.0, 5.0, 1.0])
    applylimitrule!(subj, LimitRule(lloq = 5.0, btmax = 0.0, atmax = 0.0))
    obs = MetidaNCA.getobs(subj)
    @test cmpshow("№19a C == LLOQ остаётся измеримым", 5.0, obs[3])
    @test cmpshow("№19b C < LLOQ обнуляется (контроль)", 0.0, obs[4])
end


# =============================================================================
# №20–21  ctmax
# =============================================================================
@testset "№20-21  ctmax" begin
    println("\n--- №20-21 ---------------------------------------------------------")

    # №20 при taulastp == индексу первого измеримого возвращается first(time)
    @test cmpshow("№20 tmax при ведущих NaN", 2.0,
                  MetidaNCA.tmax([0.0, 1.0, 2.0], [NaN, NaN, 7.0]))

    # №21 ctmax(::PKSubject) вычисляет taulastp, но не использует его
    subj = pkimport(pk1, :Time, :Concentration;
                    dosetime = DoseTime(dose = 100, time = 0.0, tau = 0.75))
    r  = nca!(subj)
    ct = MetidaNCA.ctmax(subj)
    @test cmpshow("№21a ctmax(subj) учитывает τ-окно", r[:Cmax], ct[1])
    @test cmpshow("№21b ctmax(subj) возвращает 4 значения, как векторный метод", 4, length(ct))
end


# =============================================================================
# №22  Отчётные границы Kel не совпадают с фактически использованными точками
# =============================================================================
@testset "№22  KelData.s/e не учитывают выброшенные нулевые точки" begin
    println("\n--- №22 ------------------------------------------------------------")

    t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    c = [0.0, 10.0, 8.0, 5.0, 0.0, 2.0]      # ноль ровно на границе диапазона (точка 5)
    subj = pkimport(t, c; kelauto = false, elimrange = ElimRange(kelstart = 3, kelend = 5))
    r    = nca!(subj)
    kd   = getkeldata(r)

    @test cmpshow("№22a число точек регрессии", 2, Int(r[:NpLZ]))
    @test cmpshow("№22b конец диапазона Kel = последняя использованная точка", 3.0, kd.e[1])
end


# =============================================================================
# №23–24  upkimport
# =============================================================================
@testset "№23-24  upkimport" begin
    println("\n--- №23-24 ---------------------------------------------------------")

    # №23a names(data) — API DataFrames, а не Tables.jl
    tbl = (subj = [1, 1, 1], conc = [1.0, 2.0, 2.0],
           st = [0.0, 1.0, 2.0], et = [1.0, 2.0, 6.0], vol = [1.0, 2.0, 3.0])
    @test nothrows("№23a upkimport из NamedTuple-таблицы",
                   () -> upkimport(tbl, :st, :et, :conc, :vol, :subj))

    # №23b проверка колонок мертва: `isa(stime, Symbol) || stime in tnames || error(...)`
    #      для Symbol всегда истинна, поэтому осмысленного сообщения нет
    got = try
        upkimport(upkdata, :no_such_column, :et, :conc, :vol, :subj); "нет ошибки"
    catch e
        _short(e)
    end
    @test cmpshow("№23b понятное сообщение об отсутствующей колонке",
                  true, occursin("Start Time", string(got)))

    # №24 promote_type(…stime…, …stime…) вместо (…stime…, …etime…)
    u24 = upkimport([0, 1, 2], [1.0, 2.0, 6.0], [1.0, 2.0, 2.0], [1.0, 2.0, 3.0])
    @test cmpshow("№24 тип времени дозы продвинут по обеим колонкам",
                  Float64, typeof(getdosetime(u24).time))
end


# =============================================================================
# №25–27, 30  Конструкторы и валидация настроек
# =============================================================================
@testset "№25-27,30  конструкторы" begin
    println("\n--- №25-27, 30 -----------------------------------------------------")

    # №25 внутренний конструктор PKSubject вызывает сам себя, если id — не Dict{Symbol,V}
    if RUN_CRASH_TESTS
        @test nothrows("№25 pkimport с id::Dict{String,Int}",
                       () -> pkimport([0.0, 1.0, 2.0], [0.0, 5.0, 2.0]; id = Dict("a" => 1)))
    else
        println("№25  пропущен (StackOverflowError). Включить: ENV[\"METIDANCA_CRASH_TESTS\"]=\"1\"")
        println("     pkimport([0.,1.,2.], [0.,5.,2.]; id = Dict(\"a\" => 1))")
        println("     Ожидается: MethodError или нормальное создание субъекта.")
        println("     Фактически: PKSubject(...,id) без типа id вызывает сам себя → StackOverflowError.")
    end

    # №26 keyword-конструктор ElimRange с вещественными границами не собирается
    #     (дефолт kelexcl = Int[] не совместим с сигнатурой (::T, ::T, ::Vector{T}))
    @test nothrows("№26 ElimRange(kelstart = 1.5, kelend = 10.0, time = true)",
                   () -> ElimRange(kelstart = 1.5, kelend = 10.0, time = true))

    # №27 setkelauto! требует kelstart > 0 даже для временного диапазона
    s27 = pkimport([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.0, 10.0, 8.0, 5.0, 3.0, 2.0])
    @test nothrows("№27 setkelrange! c ElimRange(0.0, 5.0; time = true)",
                   () -> setkelrange!(s27, ElimRange(0.0, 5.0, Float64[]; time = true); kelauto = false))

    # №30 конструкторы DoseTime вычисляют time_/dose_/tau_ и не используют их
    dt = DoseTime(dose = 100, time = 0, tau = 5)
    @test cmpshow("№30a тип time приведён к Float", Float64, typeof(dt.time))
    @test cmpshow("№30b тип dose приведён к Float", Float64, typeof(dt.dose))
    @test cmpshow("№30c тип tau  приведён к Float", Float64, typeof(dt.tau))
end


# =============================================================================
# №28, 34, 42, 45  Согласованность API
# =============================================================================
@testset "№28,34,42,45  API" begin
    println("\n--- №28, 34, 42, 45 ------------------------------------------------")

    # №28 Cltau появляется при dose = NaN, тогда как Cllast/Clinf — нет
    s28 = pkimport([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [0.0, 10.0, 8.0, 5.0, 3.0, 2.0];
                   dosetime = DoseTime(time = 0.0, tau = 4.0))     # dose = NaN
    r28 = nca!(s28)
    @test cmpshow("№28a Cllast отсутствует при dose = NaN (контроль)",
                  false, haskey(r28.result, :Cllast))
    @test cmpshow("№28b Cltau тоже должен отсутствовать",
                  false, haskey(r28.result, :Cltau))

    # №34 экспортируемый getkeldata(::PKSubject) всегда бросает ошибку
    s34 = pkimport(pk1, :Time, :Concentration)
    @test nothrows("№34 getkeldata(::PKSubject)", () -> getkeldata(s34))

    # №42 metida_table_all падает, если у субъектов разные ключи id
    a = pkimport([0.0, 1.0, 2.0], [0.0, 5.0, 2.0]; id = Dict{Symbol,Any}(:subj => 1))
    b = pkimport([0.0, 1.0, 2.0], [0.0, 4.0, 1.0]; id = Dict{Symbol,Any}(:subj => 2, :form => "T"))
    dsmix = MetidaNCA.DataSet([a, b])
    @test nothrows("№42 DataFrame(DataSet) с разными ключами id",
                   () -> DataFrame(dsmix))

    # №45 wm = true с пользовательским io (не IOBuffer)
    dswm = pkimport(pkdata2, :Time, :Concentration, :Subject)
    tmp  = tempname()
    @test nothrows("№45 nca!(ds; wm = true, io = <файл>)",
                   () -> open(f -> nca!(dswm; wm = true, io = f), tmp, "w"))
    isfile(tmp) && rm(tmp; force = true)
end


# =============================================================================
# Сводка
# =============================================================================
println("\n", "="^100)
let failed = [l for (l, ok) in SUMMARY if !ok], total = length(SUMMARY)
    printstyled("Проверок: $total, воспроизведено расхождений: $(length(failed))\n";
                color = length(failed) > 0 ? :red : :green)
    for l in failed
        println("  • ", l)
    end
end
println("="^100)
