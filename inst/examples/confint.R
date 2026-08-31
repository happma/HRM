# using the EEG dataset
z <- hrm_test(value ~ dimension, subject = "subject", data = EEG)

# simultaneous confidence intervals for each factor level combination
ci <- confint(z)
head(ci)

# which form of multiplicity control was used
attr(ci, "status")

\donttest{
# With a whole-plot factor the correlation matrix can become large enough
# that mvtnorm::qmvnorm() fails. confint() then falls back to controlling
# the family-wise error rate within each group, or to a Sidak correction.
# The status says which of these was actually used.
z <- hrm_test(value ~ group*region*variable, subject = "subject", data = EEG)
ci <- confint(z, level = 0.99)
attr(ci, "status")
}
