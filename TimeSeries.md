---
name: TimeSeries
topic: Time Series Analysis
maintainer: Rob J Hyndman, Rebecca Killick
email: Rob.Hyndman@monash.edu
version: 2022-05-04
source: https://github.com/cran-task-views/TimeSeries/
---

Base R ships with a lot of functionality useful for time series, in
particular in the stats package. This is complemented by many packages
on CRAN, which are briefly summarized below. There is overlap between
the tools for time series and those designed for specific domains including
`r view("Econometrics")`, `r view("Finance")` and `r view("Environmetrics")`.

The packages in this view can be roughly structured into the
following topics. If you think that some package is missing from the
list, please let us know, either via e-mail to the maintainer or by
submitting an issue or pull request in the GitHub repository linked above.

### Basics

-   *Infrastructure* : Base R contains substantial infrastructure for
    representing and analyzing time series data. The fundamental class
    is `"ts"` that can represent regularly spaced time series (using
    numeric time stamps). Hence, it is particularly well-suited for
    annual, monthly, quarterly data, etc.
-   *Rolling statistics* : Moving averages are computed by `ma` from
    `r pkg("forecast", priority = "core")`, and `rollmean`
    from `r pkg("zoo", priority = "core")`. The latter also
    provides a general function `rollapply`, along with other specific
    rolling statistics functions. `r pkg("slider")`
    calculates a diverse and comprehensive set of type-stable running
    functions for any R data types.
    `r pkg("tsibble", priority = "core")` provides `slide()`
    for rolling statistics, `tile()` for non-overlapping sliding
    windows, and `stretch()` for expanding windows.
    `r pkg("tbrf")` provides rolling functions based on date
    and time windows instead of n-lagged observations.
    `r pkg("roll")` provides parallel functions for
    computing rolling statistics. `r pkg("runner")` provides
    tools for running any R function in rolling windows or date windows.
    `r pkg("runstats")` provides fast computational methods
    for some running sample statistics. For
    `r pkg("data.table")`, `froll()` can be used for
    high-performance rolling statistics.
-   *Graphics* : Time series plots are obtained with `plot()` applied to
    `ts` objects. (Partial) autocorrelation functions plots are
    implemented in `acf()` and `pacf()`. Alternative versions are
    provided by `Acf()` and `Pacf()` in `r pkg("forecast")`,
    along with a combination display using `tsdisplay()`. Seasonal
    displays are obtained using `monthplot()` in stats, `seasonplot` in
    `r pkg("forecast")`, and `seasplot` in
    `r pkg("tsutils")`.
    `r pkg("feasts", priority = "core")` and
    `r pkg("brolgar")` provide various time series graphics
    for tsibble objects including time plots, season plots, subseries
    plots, ACF and PACF plots, and some combination displays.
    Interactive graphics for tsibbles using htmlwidgets are provided by
    `r pkg("tsibbletalk")`. `r pkg("SDD")`
    provides more general serial dependence diagrams, while
    `r pkg("dCovTS")` computes and plots the distance
    covariance and correlation functions of time series.
    `r pkg("ggseas")` provides additional ggplot2 graphics
    for seasonally adjusted series and rolling statistics. Calendar
    plots are implemented in `r pkg("sugrrants")`.
    `r pkg("gravitas")` allows for visualizing probability
    distributions conditional on bivariate temporal granularities.
    `r pkg("dygraphs")` provides an interface to the
    Dygraphs interactive time series charting library.
    `r pkg("TSstudio")` provides some interactive
    visualization tools for time series. `r pkg("ZRA")`
    plots forecast objects from the `r pkg("forecast")`
    package using dygraphs. Basic fan plots of forecast distributions
    are provided by `r pkg("forecast")` and
    `r pkg("vars")`. More flexible fan plots of any
    sequential distributions are implemented in
    `r pkg("fanplot")`.

### Times and Dates

-   Class `"ts"` can only deal with numeric time stamps, but many more
    classes are available for storing time/date information and
    computing with it. For an overview see *R Help Desk: Date and Time
    Classes in R* by Gabor Grothendieck and Thomas Petzoldt in [R News
    4(1)](http://CRAN.R-project.org/doc/Rnews/Rnews_2004-1.pdf), 29-32.
-   Classes `"yearmon"` and `"yearqtr"` from `r pkg("zoo")`
    allow for more convenient computation with monthly and quarterly
    observations, respectively.
-   Class `"Date"` from the base package is the basic class for dealing
    with dates in daily data. The dates are internally stored as the
    number of days since 1970-01-01.
-   The `r pkg("chron")` package provides classes for
    `dates()`, `hours()` and date/time (intra-day) in `chron()`. There
    is no support for time zones and daylight savings time. Internally,
    `"chron"` objects are (fractional) days since 1970-01-01.
-   Classes `"POSIXct"` and `"POSIXlt"` implement the POSIX standard for
    date/time (intra-day) information and also support time zones and
    daylight savings time. However, the time zone computations require
    some care and might be system-dependent. Internally, `"POSIXct"`
    objects are the number of seconds since 1970-01-01 00:00:00 GMT.
    Package `r pkg("lubridate")` provides functions that
    facilitate certain POSIX-based computations, while
    `r pkg("clock")` provides a comprehensive library for
    date-time manipulations using a new family of orthogonal date-time
    classes (durations, time points, zoned-times, and calendars).
    `r pkg("timechange")` allows for efficient manipulation
    of date-times accounting for time zones and daylight saving times.
    `r pkg("wktmo")` converts weekly data to monthly data in
    several different ways.
-   Class `"timeDate"` is provided in the
    `r pkg("timeDate")` package (previously: fCalendar). It
    is aimed at financial time/date information and deals with time
    zones and daylight savings times via a new concept of "financial
    centers". Internally, it stores all information in `"POSIXct"` and
    does all computations in GMT only. Calendar functionality, e.g.,
    including information about weekends and holidays for various stock
    exchanges, is also included.
-   The `r pkg("tis")` package provides the `"ti"` class for
    time/date information.
-   The `"mondate"` class from the `r pkg("mondate")`
    package facilitates computing with dates in terms of months.
-   The `r pkg("tempdisagg")` package includes methods for
    temporal disaggregation and interpolation of a low frequency time
    series to a higher frequency series.
-   Time series disaggregation is also provided by
    `r pkg("tsdisagg2")` and `r pkg("disaggR")`.
-   `r pkg("TimeProjection")` extracts useful time
    components of a date object, such as day of week, weekend, holiday,
    day of month, etc, and put it in a data frame.

### Time Series Classes

-   As mentioned above, `"ts"` is the basic class for regularly spaced
    time series using numeric time stamps.
-   The `r pkg("zoo")` package provides infrastructure for
    regularly and irregularly spaced time series using arbitrary classes
    for the time stamps (i.e., allowing all classes from the previous
    section). It is designed to be as consistent as possible with
    `"ts"`.
-   The package `r pkg("xts")` is based on
    `r pkg("zoo")` and provides uniform handling of R's
    different time-based data classes.
-   Several packages aim to handle time-based tibbles:
    `r pkg("tsibble")` provides tidy temporal data frames
    and associated tools; `r pkg("tsbox")` contains tools
    for working with and coercing between many time series classes
    including tsibble, ts, xts, zoo and more.
    `r pkg("timetk")` is another toolkit for converting
    between various time series data classes.
-   Some manipulation tools for time series are available in
    `r pkg("data.table")` including `shift()` for lead/lag
    operations. Further basic time series functionalities are offered by
    `r pkg("DTSg")` which is based on `r pkg("data.table")`.
    `r pkg("dtts")` provides high-frequency time series support via `r pkg("nanotime")` and `r pkg("data.table")`.
-   `r pkg("collapse")` provides fast computation of several
    time series functions such as lead/lag operations, (quasi-, log-)
    differences and growth rates on time-series and panel data, and
    ACF/PACF/CCF estimation for panel data.
-   Various packages implement irregular time series based on
    `"POSIXct"` time stamps, intended especially for financial
    applications. These include `"irts"` from
    `r pkg("tseries", priority = "core")`, and `"fts"` from
    `r pkg("fts")`.
-   The class `"timeSeries"` in `r pkg("timeSeries")`
    (previously: fSeries) implements time series with `"timeDate"` time
    stamps.
-   The class `"tis"` in `r pkg("tis")` implements time
    series with `"ti"` time stamps.
-   The package `r pkg("tframe")` contains infrastructure
    for setting time frames in different formats.
-   `r pkg("timeseriesdb")` manages time series for official
    statistics by mapping `ts` objects to PostgreSQL relations.

### Forecasting and Univariate Modeling

-   The `r pkg("fable", priority = "core")` package provides
    tools for fitting univariate time series models to many series
    simultaneously including ETS, ARIMA, TSLM and other models. It also
    provides many functions for computing and analysing forecasts. The
    time series must be in the `tsibble` format.
    `r pkg("fabletools")` provides tools for extending the
    `r pkg("fable")` framework.
-   The `r pkg("forecast")` package provides similar tools
    for `ts` objects, while `r pkg("modeltime")` and
    `r pkg("modeltime.ensemble")` provides time series
    forecasting tools for use with the 'tidymodels' ecosystem.
-   *Exponential smoothing* : `HoltWinters()` in stats provides some
    basic models with partial optimization, `ETS()` from
    `r pkg("fable")` and `ets()` from
    `r pkg("forecast")` provide a larger set of models and
    facilities with full optimization. `r pkg("robets")`
    provides a robust alternative to the `ets()` function.
    `r pkg("smooth")` implements some generalizations of
    exponential smoothing. `r pkg("legion")` implements
    multivariate versions of exponential smoothing. The
    `r pkg("MAPA")` package combines exponential smoothing
    models at different levels of temporal aggregation to improve
    forecast accuracy. Some Bayesian extensions of exponential smoothing
    are contained in `r pkg("Rlgt")`.
-   `r pkg("prophet")` forecasts time series based on an
    additive model where nonlinear trends are fit with yearly and weekly
    seasonality, plus holidays. It works best with daily data.
    `r pkg("fable.prophet")` allows prophet models to be
    used in the `r pkg("fable")` framework.
-   The theta method is implemented in the `THETA()` function from
    `r pkg("fable")`, `thetaf()` function from
    `r pkg("forecast")`, and `theta()` from
    `r pkg("tsutils")`. An alternative and extended
    implementation is provided in `r pkg("forecTheta")`.
-   *Autoregressive models* : `ar()` in stats (with model selection) and
    `r pkg("FitAR")` for subset AR models.
-   *ARIMA models* : `arima()` in stats is the basic function for ARIMA,
    SARIMA, RegARIMA, and subset ARIMA models. It is enhanced in the
    `r pkg("fable")` package via the `ARIMA()` function
    which allows for automatic modelling. Similar functionality is
    provided in the `r pkg("forecast")` package via the
    `auto.arima()` function. `arma()` in the
    `r pkg("tseries")` package provides different algorithms
    for ARMA and subset ARMA models. Other estimation methods including
    the innovations algorithm are provided by
    `r pkg("itsmr")`. Package
    `r pkg("gsarima")` contains functionality for
    Generalized SARIMA time series simulation. Robust ARIMA modeling is
    provided in the `r pkg("robustarima")` package.
    `r pkg("bayesforecast")` fits Bayesian time series
    models including seasonal ARIMA and ARIMAX models.
    `r pkg("BayesARIMAX")` implements Bayesian estimation of
    ARIMAX models. The `r pkg("mar1s")` package handles
    multiplicative AR(1) with seasonal processes.
    `r pkg("TSTutorial")` provides an interactive tutorial
    for Box-Jenkins modelling. Improved prediction intervals for ARIMA
    and structural time series models are provided by
    `r pkg("tsPI")`.
-   ARIMA models with multiple seasonal periods can be handled with
    `r pkg("tfarima")` and `r pkg("smooth")`.
-   *Periodic ARMA models* : `r pkg("partsm")` for periodic
    autoregressive time series models, and
    `r pkg("perARMA")` and `r pkg("pcts")` for
    periodic ARMA modelling and other procedures for periodic time
    series analysis.
-   *Long memory models* : Some facilities for fractional differenced
    ARFIMA models are provided in the `r pkg("fracdiff")`
    package. The `r pkg("arfima")` package has more advanced
    and general facilities for ARFIMA and ARIMA models, including
    dynamic regression (transfer function) models. Additional methods
    for fitting and simulating non-stationary ARFIMA models are in
    `r pkg("nsarfima")`. `r pkg("LongMemoryTS")`
    provides a collection of functions for analysing long memory time
    series. Fractionally differenced Gegenbaur ARMA processes are
    handled by `r pkg("garma")`.
    `r pkg("esemifar")` provides tools for nonparametric
    smoothing of long-memory time series.
-   *Transfer function* models are provided by the `arfima` function in
    the `r pkg("arfima")` and the
    `r pkg("tfarima")` packages.
-   *Structural (or unobserved component) models* are implemented in
    `StructTS()` in stats and in `r pkg("stsm")`, while
    automatic modelling and forecasting are provided by
    `r pkg("UComp")` and `r pkg("autostsm")`.
    `r pkg("KFKSDS")` provides a naive implementation of the
    Kalman filter and smoothers for univariate state space models.
    `r pkg("statespacer")` implements univariate state space
    models including structural and SARIMA models. Bayesian structural
    time series models are implemented in `r pkg("bsts")`
    Robust Kalman filtering is provided by `r pkg("RobKF")`.
-   *Non-Gaussian time series* can be handled with GLARMA state space
    models via `r pkg("glarma")`, and using Generalized
    Autoregressive Score models in the `r pkg("GAS")`
    package. `r pkg("GlarmaVarSel")` provides variable
    selection in high-dimensional sparse GLARMA models. Conditional
    Auto-Regression models using Monte Carlo Likelihood methods are
    implemented in `r pkg("mclcar")`. Efficient Bayesian
    inference for nonlinear and non-Gaussian state space models is
    provided in `r pkg("bssm")`. Non-Gaussian state space
    models with exact marginal likelihood are given by
    `r pkg("NGSSEML")`.
    `r pkg("PTSR")` includes functions to model and forecast a range of regression based dynamic models for positive time series.
-   *Count time series* models are handled in the
    `r pkg("tscount")` and `r pkg("acp")`
    packages. `r pkg("ZIM")` provides for Zero-Inflated
    Models for count time series. `r pkg("tsintermittent")`
    implements various models for analysing and forecasting intermittent
    demand time series.
-   *GARCH models* : `garch()` from `r pkg("tseries")` fits
    basic GARCH models. Many variations on GARCH models are provided by
    `r pkg("rugarch")`. Other univariate GARCH packages
    include `r pkg("fGarch")` which implements ARIMA models
    with a wide class of GARCH innovations.
    `r pkg("bayesforecast")` fits Bayesian time series
    models including several variations of GARCH models. There are many
    more GARCH packages described in the `r view("Finance")`
    task view.
-   *Stochastic volatility* models are handled by
    `r pkg("stochvol")` in a Bayesian framework.
-   *Censored time series* can be modelled using
    `r pkg("carx")`. `r pkg("ARCensReg")` fits
    univariate censored regression models with autoregressive errors.
-   *Diffusion models* such as Bass and Gompertz curves are provided by
    `r pkg("diffusion")`.
-   *Portmanteau tests* are provided via `Box.test()` in the stats
    package. Additional tests are given by
    `r pkg("portes")`,
    `r pkg("WeightedPortTest")` and
    `r pkg("testcorr")`.
-   Outlier detection following the Chen-Liu approach is provided by
    `r pkg("tsoutliers")`. The `tsoutliers` and `tsclean`
    functions in the `r pkg("forecast")` package provide
    some simple heuristic methods for identifying and correcting
    outliers. `r pkg("otsad")` implements a set of online
    anomaly detectors for time series. `r pkg("tsrobprep")`
    provides methods for replacing missing values and outliers using a
    model-based approach.
-   *Change point detection* is provided in
    `r pkg("strucchange")` and
    `r pkg("strucchangeRcpp")` (using linear regression
    models) and in `r pkg("trend")` (using nonparametric
    tests). The `r pkg("changepoint")` package provides many
    popular changepoint methods, and `r pkg("ecp")` does
    nonparametric changepoint detection for univariate and multivariate
    series. `r pkg("changepoint.np")` implements the
    nonparametric PELT algorithm, `r pkg("changepoint.mv")`
    detects changepoints in multivariate time series, while
    `r pkg("changepoint.geo")` implements the
    high-dimensional changepoint detection method GeomCP.
    Factor-augmented VAR (FAVAR) models are estimated by a Bayesian method with `r pkg("FAVAR")`.
    `r pkg("InspectChangepoint")` uses sparse projection to
    estimate changepoints in high-dimensional time series.
    `r pkg("Rbeast")` provides Bayesian change-point detection and time series decomposition.
    `r pkg("breakfast")` includes methods for fast multiple
    change-point detection and estimation.
-   Tests for possibly non-monotonic trends are provided by
    `r pkg("funtimes")`.
-   *Time series imputation* is provided by the
    `r pkg("imputeTS")` package. Some more limited
    facilities are available using `na.interp()` from the
    `r pkg("forecast")` package.
    `r pkg("imputeTestbench")` provides tools for testing
    and comparing imputation methods. `r pkg("mtsdi")`
    implements an EM algorithm for imputing missing values in
    multivariate normal time series, accounting for spatial and temporal
    correlations.
-   The `r pkg("seer")` package implements a framework for
    feature-based forecast model selection.
-   Forecasts can be combined in the `r pkg("fable")`
    package using simple linear expressions.
    `r pkg("ForecastComb")` supports many forecast
    combination methods including simple, geometric and regression-based
    combinations. `r pkg("forecastHybrid")` provides
    functions for ensemble forecasts, combining approaches from the
    `r pkg("forecast")` package.
    `r pkg("opera")` has facilities for online predictions
    based on combinations of forecasts provided by the user.
    `r pkg("profoc")` combines probabilistic forecasts using
    CRPS learning.
-   Point forecast evaluation is provided in the `accuracy()` function
    from the `r pkg("fable")` and
    `r pkg("forecast")` packages. Distributional forecast
    evaluation using scoring rules is available in
    `r pkg("fable")`, `r pkg("scoringRules")`
    and `r pkg("scoringutils")`. The Diebold-Mariano test
    for comparing the forecast accuracy of two models is implemented in
    the `dm.test()` function in `r pkg("forecast")`.
    `r pkg("tsutils")` implements the Nemenyi test for comparing forecasts.
    `r pkg("greybox")` provides `ro()` for general rolling
    origin evaluation of forecasts.
-   Tidy tools for forecasting are provided by
    `r pkg("sweep")`, converting objects produced in
    `r pkg("forecast")` to "tidy" data frames.
-   Multi-step-ahead direct forecasting with several machine learning
    approaches are provided in `r pkg("forecastML")`.
-   *Miscellaneous* : `r pkg("ltsa")` contains methods for
    linear time series analysis, `r pkg("timsac")` for time
    series analysis and control.

### Frequency analysis

-   *Spectral density estimation* is provided by `spectrum()` in the
    stats package, including the periodogram, smoothed periodogram and
    AR estimates. Bayesian spectral inference is provided by
    `r pkg("bspec")`, `r pkg("beyondWhittle")` and `r pkg("regspec")`.
    `r pkg("quantspec")` includes methods to compute and
    plot Laplace periodograms for univariate time series. The
    Lomb-Scargle periodogram for unevenly sampled time series is
    computed by `r pkg("lomb")`.
    `r pkg("spectral")` uses Fourier and Hilbert transforms
    for spectral filtering. `r pkg("psd")` produces
    adaptive, sine-multitaper spectral density estimates.
    `r pkg("kza")` provides Kolmogorov-Zurbenko Adaptive
    Filters including break detection, spectral analysis, wavelets and
    KZ Fourier Transforms. `r pkg("multitaper")` also
    provides some multitaper spectral analysis tools. Higher-order
    spectral analysis is implemented in `r pkg("rhosa")`,
    including bispectrum, bicoherence, cross-bispectrum and
    cross-bicoherence.
-   *Wavelet methods* : The `r pkg("wavelets")` package
    includes computing wavelet filters, wavelet transforms and
    multiresolution analyses. Multiresolution forecasting using wavelets
    is also implemented in `r pkg("mrf")`.
    `r pkg("WaveletComp")` provides some tools for
    wavelet-based analysis of univariate and bivariate time series
    including cross-wavelets, phase-difference and significance tests.
    `r pkg("biwavelet")` is a port of the WTC Matlab package for univariate and bivariate wavelet analyses.
    `r pkg("mvLSW")` provides tools for multivariate locally stationary wavelet processes.
    `r pkg("LSWPlib")` contains functions for simulation and spectral estimation of locally stationary wavelet packet processes.
    Tests of white noise using wavelets are provided by `r pkg("hwwntest")`.
    Wavelet scalogram tools are contained in `r pkg("wavScalogram")`.
    Further wavelet methods can be found in the packages
    `r pkg("rwt")`, `r pkg("waveslim")`,
    `r pkg("wavethresh")`.
-   *Harmonic regression* using Fourier terms is implemented in
    `r pkg("HarmonicRegression")`. The
    `r pkg("fable")` and `r pkg("forecast")`
    packages also provide some simple harmonic regression facilities via
    the `fourier` function.

### Decomposition and Filtering

-   *Filters and smoothing* : `filter()` in stats provides
    autoregressive and moving average linear filtering of multiple
    univariate time series. The `r pkg("robfilter")` package
    provides several robust time series filters. `smooth()` from the
    stats package computes Tukey's running median smoothers, 3RS3R,
    3RSS, 3R, etc. `r pkg("sleekts")` computes the 4253H
    twice smoothing method. `r pkg("mFilter")` implements
    several filters for smoothing and extracting trend and cyclical
    components including Hodrick-Prescott and Butterworth filters.
    `r pkg("smoots")` provides nonparametric estimation of
    the time trend and its derivatives.
-   *Decomposition* : Seasonal decomposition is discussed below.
    Autoregressive-based decomposition is provided by
    `r pkg("ArDec")`. `r pkg("tsdecomp")`
    implements ARIMA-based decomposition of quarterly and monthly data.
    `r pkg("rmaf")` uses a refined moving average filter for
    decomposition.
-   *Singular Spectrum Analysis* is implemented in
    `r pkg("Rssa")` and `r pkg("ASSA")`.
-   *Empirical Mode Decomposition* (EMD) and Hilbert spectral analysis
    is provided by `r pkg("EMD")`. Additional tools,
    including ensemble EMD, are available in `r pkg("hht")`.
    An alternative implementation of ensemble EMD and its complete
    variant are available in `r pkg("Rlibeemd")`.

### Seasonality

-   *Seasonal decomposition* : the stats package provides classical
    decomposition in `decompose()`, and STL decomposition in `stl()`.
    Enhanced STL decomposition is available in
    `r pkg("stlplus")`. `r pkg("stR")` provides
    Seasonal-Trend decomposition based on Regression.
    `r pkg("smooth")` and `r pkg("tsutils")`
    implement extended versions of classical decomposition.
-   X-13-ARIMA-SEATS binaries are provided in the
    `r pkg("x13binary")` package, with
    `r pkg("seasonal")` providing an R interface and
    `r pkg("seasonalview")` providing a GUI. An alternative
    interface is provided by `r pkg("x12")`.
-   An interface to the JDemetra+ seasonal adjustment software is
    provided by `r pkg("RJDemetra")`.
    `r pkg("ggdemetra")` provides associated ggplot2
    functions.
-   Seasonal adjustment of daily time series, allowing for day-of-week,
    time-of-month, time-of-year and holiday effects is provided by
    `r pkg("dsa")`.
-   *Analysis of seasonality* : the `r pkg("bfast")` package
    provides methods for detecting and characterizing abrupt changes
    within the trend and seasonal components obtained from a
    decomposition. `r pkg("npst")` provides a generalization
    of Hewitt's seasonality test.
-   `r pkg("season")`: Seasonal analysis of health data
    including regression models, time-stratified case-crossover,
    plotting functions and residual checks.
-   `r pkg("seas")`: Seasonal analysis and graphics,
    especially for climatology.
-   `r pkg("deseasonalize")`: Optimal deseasonalization for
    geophysical time series using AR fitting.
-   `r pkg("sazedR")`: Method to estimate the period of a
    seasonal time series.

### Stationarity, Unit Roots, and Cointegration

-   *Stationarity and unit roots* : `r pkg("tseries")`
    provides various stationarity and unit root tests including
    Augmented Dickey-Fuller, Phillips-Perron, and KPSS. Alternative
    implementations of the ADF and KPSS tests are in the
    `r pkg("urca")` package, which also includes further
    methods such as Elliott-Rothenberg-Stock, Schmidt-Phillips and
    Zivot-Andrews tests. `r pkg("uroot")` provides seasonal
    unit root tests. `r pkg("CADFtest")` provides
    implementations of both the standard ADF and a covariate-augmented
    ADF (CADF) test. `r pkg("MultipleBubbles")` tests for
    the existence of bubbles based on Phillips-Shi-Yu (2015).
-   *Local stationarity* : `r pkg("locits")` provides a test
    of local stationarity and computes the localized autocovariance.
    Time series costationarity determination is provided by
    `r pkg("costat")`. `r pkg("LSTS")` has
    functions for locally stationary time series analysis. Locally
    stationary wavelet models for nonstationary time series are
    implemented in `r pkg("wavethresh")` (including
    estimation, plotting, and simulation functionality for time-varying
    spectra).
-   *Cointegration* : The Engle-Granger two-step method with the
    Phillips-Ouliaris cointegration test is implemented in
    `r pkg("tseries")` and `r pkg("urca")`. The
    latter additionally contains functionality for the Johansen trace
    and lambda-max tests. `r pkg("tsDyn")` provides
    Johansen's test and AIC/BIC simultaneous rank-lag selection.
    Parameter estimation and inference in a cointegrating regression are
    implemented in `r pkg("cointReg")`.
    `r pkg("nardl")` estimates nonlinear cointegrating
    autoregressive distributed lag models.

### Nonlinear Time Series Analysis

-   *Nonlinear autoregression* : Tools for nonlinear time series
    analysis are provided in `r pkg("NTS")` including
    threshold autoregressive models, Markov-switching models,
    convolutional functional autoregressive models, and nonlinearity
    tests. Various forms of nonlinear autoregression are available in
    `r pkg("tsDyn")` including additive AR, SETAR and LSTAR models,
    threshold VAR and VECM.
    `r pkg("bentcableAR")` implements Bent-Cable autoregression.
    `r pkg("BAYSTAR")` provides Bayesian analysis of threshold autoregressive models.
    Mixture AR models are implemented in `r pkg("mixAR")` and `r pkg("uGMAR")`.
-   *Neural network autoregression* : Neural network forecasting based on
    lagged inputs are provided by `r pkg("tsDyn")`, `r pkg("GMDH")` and
    `r pkg("nnfor")`. `r pkg("NlinTS")` includes neural
    network VAR, and a nonlinear version of the Granger causality test
    based on feedforward neural networks.
    `r pkg("TSLSTM")` provides forecasts using a Long Short Term Memory (LSTM) model.
    `r pkg("TSANN")` automatically identifies an artificial neural network
    based on forecasting accuracy.
-   `r pkg("tseriesChaos")` provides an R implementation of
    the algorithms from the
    *[TISEAN](http://www.mpipks-dresden.mpg.de/~tisean/) project* .
    `r pkg("DChaos")` provides several algorithms for
    detecting chaotic signals inside univariate time series.
-   Autoregression Markov switching models are provided in
    `r pkg("MSwM")`, while dependent mixtures of latent
    Markov models are given in `r pkg("depmix")` and
    `r pkg("depmixS4")` for categorical and continuous time
    series.
-   *Tests* : Various tests for nonlinearity are provided in
    `r pkg("fNonlinear")`.
    `r pkg("tseriesEntropy")` tests for nonlinear serial
    dependence based on entropy metrics.
-   Additional functions for nonlinear time series are available in
    `r pkg("nlts")` and
    `r pkg("nonlinearTseries")`.

### Entropy

-   `r pkg("RTransferEntropy")` measures information flow
    between time series with Shannon and Renyi transfer entropy.
-   An entropy measure based on the Bhattacharya-Hellinger-Matusita
    distance is implemented in `r pkg("tseriesEntropy")`.
-   Various approximate and sample entropies are computed using
    `r pkg("TSEntropies")`.

### Dynamic Regression Models

-   *Dynamic linear models* : A convenient interface for fitting dynamic
    regression models via OLS is available in
    `r pkg("dynlm")`; an enhanced approach that also works
    with other regression functions and more time series classes is
    implemented in `r pkg("dyn")`. More advanced dynamic
    system equations can be fitted using `r pkg("dse")`.
    Gaussian linear state space models can be fitted using
    `r pkg("dlm")` (via maximum likelihood, Kalman
    filtering/smoothing and Bayesian methods), or using
    `r pkg("bsts")` which uses MCMC.
    `r pkg("dLagM")` provides time series regression with
    distributed lags. Functions for distributed lag nonlinear modelling
    are provided in `r pkg("dlnm")`.
    `r pkg("sym.arma")` will fit ARMA models with regressors
    where the observations follow a conditional symmetric distribution.
-   *Time-varying parameter* models can be fitted using the
    `r pkg("tpr")` package.
-   `r pkg("greybox")` provides several tools for modelling
    and forecasting with dynamic regression models.

### Multivariate Time Series Models

-   *Vector autoregressive (VAR) models* are provided via `ar()` in the
    basic stats package including order selection via the AIC. These
    models are restricted to be stationary.
     `r pkg("MTS")` is an all-purpose toolkit for analyzing multivariate time series including VAR, VARMA, seasonal VARMA, VAR models with exogenous variables, multivariate regression with time series errors, and much more.
    Possibly non-stationary VAR models are fitted in the
    `r pkg("mAr")` package, which also allows VAR models in
    principal component space. `r pkg("sparsevar")` allows
    estimation of sparse VAR and VECM models,
    `r pkg("bigtime")` estimates large sparse VAR, VARX and
    VARMA models, while `r pkg("BigVAR")` estimates VAR and
    VARX models with structured lasso penalties and
    `r pkg("svars")` implements data-driven structural VARs.
    Shrinkage estimation methods for VARs are implemented in
    `r pkg("VARshrink")`. More elaborate models are provided
    in package `r pkg("vars")`,
    `r pkg("tsDyn")`, `estVARXls()` in
    `r pkg("dse")`. Another implementation with bootstrapped
    prediction intervals is given in `r pkg("VAR.etp")`.
    `r pkg("bvartools")` assists in the set-up of Bayesian
    VAR models, while `r pkg("BMTAR")` implements Baysian
    Multivariate Threshold AR models with missing data.
    `r pkg("mfbvar")` includes tools for estimating
    mixed-frequency Bayesian VAR models. `r pkg("BVAR")`
    provides a toolkit for hierarchical Bayesian VAR models.
    `r pkg("BGVAR")` implements Bayesian Global VAR models.
    `r pkg("mlVAR")` provides multi-level vector
    autoregression. `r pkg("VARsignR")` provides routines
    for identifying structural shocks in VAR models using sign
    restrictions. `r pkg("gmvarkit")` estimates Gaussian
    mixture VAR models. `r pkg("GNAR")` provides methods for
    fitting network AR models, while `r pkg("graphicalVAR")`
    estimates graphical VAR models. `r pkg("gdpc")`
    implements generalized dynamic principal components.
    `r pkg("pcdpca")` extends dynamic principal components
    to periodically correlated multivariate time series.
    `r pkg("onlineVAR")` implements online fitting of
    time-adaptive lasso VARs. `r pkg("mgm")` estimates
    time-varying mixed graphical models and mixed VAR models via
    regularized regression.
-   *VARIMA models* and *state space models* are provided in the
    `r pkg("dse")` package.
-   *Vector error correction models* are available via the
    `r pkg("urca")`, `r pkg("ecm")`,
    `r pkg("vars")`, `r pkg("tsDyn")` packages,
    including versions with structural constraints and thresholding.
-   *Vector exponential smoothing* is provided by
    `r pkg("smooth")`.
-   *Time series component analysis* : `r pkg("ForeCA")`
    implements forecastable component analysis by searching for the best
    linear transformations that make a multivariate time series as
    forecastable as possible. `r pkg("PCA4TS")` finds a
    linear transformation of a multivariate time series giving
    lower-dimensional subseries that are uncorrelated with each other.
    `r pkg("HDTSA")` provides procedures for several
    high-dimensional time series analysis tools. One-sided dynamic
    principal components are computed in `r pkg("odpc")`.
    Frequency-domain-based dynamic PCA is implemented in
    `r pkg("freqdom")`. `r pkg("tsBSS")`
    provides blind source separation and supervised dimension reduction
    for time series.
-   *Multivariate state space models* An implementation is provided by
    the `r pkg("KFAS")` package which provides a fast
    multivariate Kalman filter, smoother, simulation smoother and
    forecasting. `r pkg("FKF")` provides a fast and flexible
    implementation of the Kalman filter, which can deal with missing
    values. `r pkg("FKF.SP")` implements fast Kalman
    filtering through sequential processing. Another implementation is
    given in the `r pkg("dlm")` package which also contains
    tools for converting other multivariate models into state space
    form. `r pkg("mssm")` also provides methods for
    multivariate state space models. `r pkg("MARSS")` fits
    constrained and unconstrained multivariate autoregressive
    state-space models using an EM algorithm.
    `r pkg("mbsts")` provides tools for multivariate
    Bayesian structural time series models. All of these packages assume
    the observational and state error terms are uncorrelated.
-   *Partially-observed Markov processes* are a generalization of the
    usual linear multivariate state space models, allowing non-Gaussian
    and nonlinear models. These are implemented in the
    `r pkg("pomp")` package.
-   Multivariate stochastic volatility models (using latent factors) are
    provided by `r pkg("factorstochvol")`.

### Analysis of large groups of time series

-   *Time series features* are computed in `r pkg("feasts")`
    for time series in `tsibble` format. They are computed using
    `r pkg("tsfeatures")` for a list or matrix of time
    series in `ts` format. In both packages, many built-in feature
    functions are included, and users can add their own.
    `r pkg("Rcatch22")` provides fast computation of 22
    features identified as particularly useful.
    `r pkg("theft")` calculates time series features from various R and Python packages.
    `r pkg("fsMTS")` implements feature selection routines
    for multivariate time series.
-   *Time series clustering* is implemented in
    `r pkg("TSclust")`, `r pkg("dtwclust")`,
    `r pkg("BNPTSclust")` and `r pkg("pdc")`.
-   `r pkg("TSdist")` provides distance measures for time
    series data.
-   `r pkg("TSrepr")` includes methods for representing time
    series using dimension reduction and feature extraction.
-   `r pkg("rucrdtw")` provides R bindings for functions
    from the UCR Suite to enable ultrafast subsequence search for a best
    match under Dynamic Time Warping and Euclidean Distance.
    `r pkg("IncDTW")` provides incremental calculation of
    dynamic time warping for streaming time series.
-   Methods for plotting and forecasting collections of hierarchical and
    grouped time series are provided by `r pkg("fable")` and
    `r pkg("hts")`. `r pkg("thief")` uses
    hierarchical methods to reconcile forecasts of temporally aggregated
    time series. `r pkg("FoReco")` provides various forecast
    reconciliation methods for cross-sectional, temporal, and
    cross-temporal constrained time series. An alternative approach to
    reconciling forecasts of hierarchical time series is provided by
    `r pkg("gtop")`. `r pkg("ProbReco")`
    provides tools to train forecast reconciliation weights by
    optimizing probability scoring functions.

### Functional time series

-   Tools for visualizing, modeling, forecasting and analysing
    functional time series are implemented in
    `r pkg("ftsa")`. `r pkg("NTS")` also
    implements functional autoregressive models. Seasonal functional
    autoregression models are provided by `r pkg("Rsfar")`.
    `r pkg("fpcb")` implements predictive confidence bands
    for functional time series.
-   `r pkg("fdaACF")` estimates the autocorrelation function
    for functional time series.
-   `r pkg("freqdom.fda")` provides implements of dynamical
    functional principal components for functional time series.
-   `r pkg("STFTS")` contains stationarity, trend and unit
    root tests for functional time series.
-   `r pkg("wwntests")` provides an array of white noise
    hypothesis tests for functional data.

### Matrix and tensor-valued time series

-   `r pkg("tensorTS")` provides functions for estimation,
    simulation and prediction of factor and autoregressive models for
    matrix and tensor valued time series.

### Continuous time models

-   `r pkg("carfima")` allows for continuous-time ARFIMA
    models.
-   `r pkg("Sim.DiffProc")` simulates and models stochastic
    differential equations.
-   Simulation and inference for stochastic differential equations is
    provided by `r pkg("sde")` and
    `r pkg("yuima")`.

### Resampling

-   *Bootstrapping* : The `r pkg("boot")` package provides
    function `tsboot()` for time series bootstrapping, including block
    bootstrap with several variants. `r pkg("blocklength")`
    allows for selecting the optimal block-length for a dependent
    bootstrap. `tsbootstrap()` from `r pkg("tseries")`
    provides fast stationary and block bootstrapping. Maximum entropy
    bootstrap for time series is available in
    `r pkg("meboot")`. `r pkg("timesboot")`
    computes the bootstrap CI for the sample ACF and periodogram.
    `r pkg("BootPR")` computes bias-corrected forecasting
    and bootstrap prediction intervals for autoregressive time series.
    `r pkg("bootUR")` implements bootstrap unit root tests.

### Time Series Data

-   Various data sets in `r pkg("tsibble")` format are
    provided by `r pkg("tsibbledata")`.
-   Data from Cryer and Chan (2010, 2nd ed) *Time series analysis with
    applications in R* are in the `r pkg("TSA")` package.
-   Data from Hyndman and Athanasopoulos (2018, 2nd ed) *Forecasting:
    principles and practice* are in the `r pkg("fpp2")`
    package.
-   Data from Hyndman and Athanasopoulos (2021, 3rd ed) *Forecasting:
    principles and practice* are in the `r pkg("fpp3")`
    package.
-   Data from Hyndman, Koehler, Ord and Snyder (2008) *Forecasting with
    exponential smoothing* are in the `r pkg("expsmooth")`
    package.
-   Data from Makridakis, Wheelwright and Hyndman (1998, 3rd ed)
    *Forecasting: methods and applications* are in the
    `r pkg("fma")` package.
-   Data from Shumway and Stoffer (2017, 4th ed) *Time Series Analysis
    and Its Applications: With R Examples* are in the
    `r pkg("astsa")` package.
-   Data from Tsay (2005, 2nd ed) *Analysis of Financial Time Series*
    are in the `r pkg("FinTS")` package.
-   Data from Woodward, Gray, and Elliott (2016, 2nd ed) *Applied Time
    Series Analysis with R* are in the `r pkg("tswge")`
    package.
-   `r pkg("AER")` and `r pkg("Ecdat")` both
    contain many data sets (including time series data) from many
    econometrics text books
-   Data from the M-competition and M3-competition are provided in the
    `r pkg("Mcomp")` package. `r pkg("Tcomp")`
    provides data from the 2010 IJF Tourism Forecasting Competition.
-   *National time series data:* `r pkg("readabs")`
    downloads, imports and tidies time series data from the
    [*Australian* Bureau of Statistics](https://www.abs.gov.au).
    `r pkg("BETS")` provides access to the most important
    economic time series in *Brazil* . `r pkg("bundesbank")`
    allows access to the time series databases of the *German* central
    bank. Data from *Switzerland* via
    [dataseries.org](http://dataseries.org) can be downloaded and
    imported using `r pkg("dataseries")`.
    `r pkg("ugatsdb")` provides an API to access time series
    data for *Uganda* .
-   *Time series data bases:* `r pkg("fame")` provides an
    interface for FAME time series databases. Economic time series and
    other data from FRED (the Federal Reserve Economic Data) can be
    retrieved using `r pkg("fredr")`.
    `r pkg("rdbnomics")` provides access to hundreds of
    millions of time series from [DBnomics](db.nomics.world).
    `r pkg("influxdbr")` provides an interface to the
    InfluxDB time series database. `r pkg("pdfetch")`
    provides facilities for downloading economic and financial time
    series from public sources. Data from the
    [Quandl](http://www.quandl.com) online portal to financial,
    economical and social datasets can be queried interactively using
    the `r pkg("Quandl")` package.
    `r pkg("TSdbi")` provides a common interface to time
    series databases. `r pkg("tsdb")` implements a simple
    data base for numerical time series.
-   *Synthetic data* are produced by `simulate()` in
    `r pkg("forecast")` package or `generate()` in
    `r pkg("fable")`, given a specific model.
    `r pkg("gratis")` generates new time series with diverse
    and controllable characteristics using mixture autoregression
    models. `r pkg("synthesis")` generates synthetic time
    series from commonly used statistical models, including linear,
    nonlinear and chaotic systems. `r pkg("tssim")` flexibly
    simulates daily or monthly time series using seasonal, calendar, and
    outlier components.

### Miscellaneous

-   `r pkg("dtw")`: Dynamic time warping algorithms for
    computing and plotting pairwise alignments between time series.
-   `r pkg("EBMAforecast")`: Ensemble Bayesian model
    averaging forecasts using Gibbs sampling or EM algorithms.
-   `r pkg("ensembleBMA")`: Bayesian Model Averaging to
    create probabilistic forecasts from ensemble forecasts and weather
    observations.
-   `r pkg("earlywarnings")`: Early warnings signals toolbox
    for detecting critical transitions in time series
-   `r pkg("FeedbackTS")`: Analysis of fragmented time
    directionality to investigate feedback in time series.
-   `r pkg("gsignal")` is an R implementation of the Octave
    package "signal", containing a variety of signal processing tools.
-   `r pkg("LPStimeSeries")` aims to find "learned pattern
    similarity" for time series.
-   `r pkg("nets")`: routines for the estimation of sparse
    long run partial correlation networks for time series data.
-   `r pkg("paleoTS")`: Modeling evolution in
    paleontological time series.
-   `r pkg("pastecs")`: Regulation, decomposition and
    analysis of space-time series.
-   `r pkg("PSF")`: Forecasting univariate time series using
    pattern-sequences.
-   `r pkg("ptw")`: Parametric time warping.
-   `r pkg("RGENERATE")` provides tools to generate vector
    time series.
-   `r pkg("RMAWGEN")` is set of S3 and S4 functions for
    spatial multi-site stochastic generation of daily time-series of
    temperature and precipitation making use of VAR models. The package
    can be used in climatology and statistical hydrology.
-   `r pkg("RSEIS")`: Seismic time series analysis tools.
-   `r pkg("rts")`: Raster time series analysis (e.g., time
    series of satellite images).
-   `r pkg("SLBDD")`: Functions for analyzing large-scale time series, based on the book "Statistical Learning with Big Dependent Data" (Pena & Tsay, 2021).
-   `r pkg("spTimer")`: Spatio-temporal Bayesian modelling.
-   `r pkg("surveillance")`: Temporal and spatio-temporal
    modeling and monitoring of epidemic phenomena.
-   `r pkg("Tides")`: Functions to calculate characteristics
    of quasi periodic time series, e.g. observed estuarine water levels.
-   `r pkg("tiger")`: Temporally resolved groups of typical
    differences (errors) between two time series are determined and
    visualized.
-   `r pkg("tsfknn")`: Time series forecasting with
    k-nearest-neighbours.
-   `r pkg("TSMining")`: Mining Univariate and Multivariate
    Motifs in Time-Series Data.
-   `r pkg("tsModel")`: Time series modeling for air
    pollution and health.

### Links
-   [TISEAN Project](http://www.mpipks-dresden.mpg.de/~tisean/)
