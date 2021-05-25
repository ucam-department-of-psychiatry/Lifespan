===== README =====

=== Scripts ===

== 1xx ===

100.common-variables.r
101.common-functions.r
102.gamlss-recode.r

Common functions and a re-write of several gamlss functions (that had
issues within our pipeline).

== 2xx ==

200.variables.r
201.functions.r
202.functions-plotting.r
210.data-setup.r
220.simulation-omega-setup.r

Import and clean the data

== 3xx ==

300.variables.r
301.functions.r
310.fitting.r
320.best-fit.r
330.bootstrapping.r
340.bootstrap-merge.r
350.calc-derived.r
350.calc-novel.r

Main scripts, these fit the gamlss model(s), select the best (via
BIC), perform bootstrapping, and calculate all necessary derived values.

== 4xx ==

400.variables.r
420.simulation-omega-copy-fit-boot-best.r

Only for the simulated datasets, used to illustrate the impact of:

FIT(SUBSET+NOVEL) - EXPAND(FIT(SUBSET),NOVEL)

== 5xx ==

500.plotting-variables.r
501.plotting-functions.r
510.plotting.r

Plotting functions, these //only// use the DERIVED.rds and the NOVEL
fitted objects.

