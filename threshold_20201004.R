# US Health Opportunity Cost Threshold Analysis
# Copyright (C) 2020 David J. Vanness
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

rm(list = ls())

vstr = "20201004" #Version
M = 50000 #Sensitivity analysis replications / set to 1 for base case results, 50000 for PSA / Note: run base case first, then PSA

# GAM estimate of SF-6D by Age -----uncomment this section to obtain GAM results for SF-6D by age / -------------------------------------------------------------------------------------
# Note: results in published paper used GAM 1.16.1 in R 3.6.1 : later versions of GAM and R may change bootstrap results trivially
#
# library(readr)
# require(haven)
# require(gam)
# require(modelr)
# require(dplyr)
#
# set.seed(seed = 112358)
# data.in <- read_dta(file = "data/DS0001/23263-0001-Data.dta") #Note: to obtain this file, you must download from: https://www.disc.wisc.edu/archivereport/downloadform.asp?studyID=70
# data.in <- select(data.in,c(AGE,SF6D_12V2,HUI3SCORE,EQ5DSCORE,NOINSURE)) %>%
#      filter(NOINSURE==2)
#
# df_ = 4
#
# model.gam <- gam(SF6D_12V2 ~ ns(x = AGE,df=df_,intercept = TRUE,Boundary.knots = c(0,100)),data = data.in)
# data.predict <- data.frame(AGE=0:100)
# predict.out <- predict(object = model.gam,newdata = data.predict)
# null_ <- write.table(x = cbind(0:100,predict.out),file = "SF6D_12V2.txt",row.names = FALSE,col.names = FALSE)
#
# boot.predict.out <- matrix(data = NA,nrow = 50000,ncol = 101)
# for (j in 1:50000){
#      boot.model.gam <- gam(SF6D_12V2 ~ ns(x = AGE,df=df_,intercept = TRUE,Boundary.knots = c(0,100)),data = resample_bootstrap(data.in))
#      boot.predict.out[j,] <- predict(object = boot.model.gam,newdata = data.predict)
# }
#
# boot.result <- apply(boot.predict.out, MARGIN = 2,quantile, probs=c(.025,.975))
# boot.result[2,] <- pmin(boot.result[2,],1)
#
# save(list = c("boot.result","boot.predict.out"),file = paste("data/sf6D_boot_",vstr,".Rdata",sep = ""))

#Declare Functions -----------------------------------------------------------------------------------------------------
findgp = function(sh, mu, lt, p) {
        pgamma(q = lt,
               shape = sh,
               scale = mu / sh) - p
} #Finds shape of gamma with mean mu and odds p/(1-p) less than cutoff lt
findbp = function(phi, mu, lt, p) {
        pbeta(lt, phi * mu, phi * (1 - mu)) - p
} #Finds precision of beta with mean mu and odds p/(1-p) less than cutoff lt
findsbp = function(phi, mu, lt, lb, ub, p) {
        #Finds precision of scaled beta with bounds [lb,ub], mean mu and odds p/(1-p) less than cutoff lt
        pbeta(
                q = (lt - lb) / (ub - lb),
                shape1 = phi * ((mu - lb) / (ub - lb)),
                shape2 = phi * ((ub - mu) / (ub - lb))
        ) - p
}

#Load data--------------------------------------------------------------------------------------------------------------
indat = read.csv(file = paste("data/popdata_", vstr, ".csv", sep = ""),
                 header = TRUE)
attach(indat)
load(file = paste("data/sf6D_boot_", vstr, ".Rdata", sep = "")) #

#Hyperparameters for Model Parameters

#P: Average Direct-Purchase Private Insurance Premiums------------------------------------------------------------------
mu_P = 12 * ((8411614 * 612) - (7325211 * (539 + 87))) / (8411164 - 7325211) #Mean exchange non-APTC premium for HC.gov states (2019 OEP State PUF, line 53)
shape_P = uniroot(
        f = findgp,
        interval = c(.01, 200),
        mu = mu_P,
        lt = .8 * mu_P,
        p = 1 / 101
)$root #Odds that the average premium is less than 80% of its base-case value are 1:100
scale_P = mu_P / shape_P

#NNT: Number needed to lose insurance to result in one statistical death------------------------------------------------
mu_NNT = 277.5 #Expected NNT from Sommers 2017 [range: 239-316]
#mu_NNT = 310 #Expected NNT from Borgschulte and Vogler 2019 - Uncomment for scenario analysis
shape_NNT = uniroot(
        f = findgp,
        interval = c(.01, 200),
        mu = mu_NNT,
        lt = .5 * mu_NNT,
        p = 1 / 101
)$root #Odds that NNT is less than 50% of its base-case value are 1:100
scale_NNT = mu_NNT / shape_NNT

#alpha: Percentage of age-related morbidity amenable to health care-----------------------------------------------------
mu_alpha = 0.10 #Kaplan and Milstein 2019
phi_alpha = uniroot(
        f = findbp,
        mu = mu_alpha,
        lt = .5 * mu_alpha,
        p = 1 / 101,
        interval = c(.01, 200)
)$root #Odds that NNT is less than half of its base-case value are 1:100

#theta: Proportion of premium passed through to policyholders-----------------------------------------------------------
mu_theta = 1 #Derived from profit maximization given monopolistic competition: MR=MC
lb_theta = .5
ub_theta = 1.5
phi_theta = uniroot(
        f = findsbp,
        interval = c(.01, 200),
        mu = mu_theta,
        lb = lb_theta,
        ub = ub_theta,
        lt = .8 * mu_theta,
        p = 1 / 101
)$root #Theta in [0.5,1.5], odds that theta is less than 80% its base-case value are 1:100

#eta: Coverage elasticity with respect to premium-----------------------------------------------------------------------
mu_etaP_j = c(-1.5,-1.05,-.7) #Age 18-34, 35-54, 55-64 from Saltzman 2019
ub_eta = 3
lb_eta = 0
phi_eta1 = uniroot(
        f = findsbp,
        interval = c(.01, 200),
        mu = -mu_etaP_j[1],
        lb = lb_eta,
        ub = ub_eta,
        lt = (1 / 3) * abs(mu_etaP_j[1]),
        p = 1 / 101
)$root #Age 18-34 - Odds that |elasticity| is less than 1/3 its base case value are 1:100
phi_eta2 = uniroot(
        f = findsbp,
        interval = c(.01, 200),
        mu = -mu_etaP_j[2],
        lb = lb_eta,
        ub = ub_eta,
        lt = (1 / 3) * abs(mu_etaP_j[2]),
        p = 1 / 101
)$root #Age 35-54 - Odds that |elasticity| is less than 1/3 its base case value are 1:100
phi_eta3 = uniroot(
        f = findsbp,
        interval = c(.01, 200),
        mu = -mu_etaP_j[3],
        lb = lb_eta,
        ub = ub_eta,
        lt = (1 / 3) * abs(mu_etaP_j[3]),
        p = 1 / 101
)$root #Age 55-64 - Odds that |elasticity| is less than 1/3 its base case value are 1:100

# Age distribution of the direct purchase private insurance population---------------------------------------------------
mu_dirpi17 = 1000 * c(2358, 2596, 3163, 4327, 5967, 5138, 6675, 7401) #2017 Estimated population covered by direct
# purchase insurance ages 0-5, 6-11, 12-17, 18-24, 25-34, 35-44, 45-54, 55-64 from Table HI01, 2018 CPS ASEC Supplement
# https://www.census.gov/data/tables/time-series/demo/income-poverty/cps-hi/hi-01.html
se_dirpi17 = sqrt(-0.000010 * mu_dirpi17 ^ 2 + 3240 * mu_dirpi17) #Generalized Variance Parameters from
# https://www2.census.gov/programs-surveys/cps/techdocs/cpsmar18.pdf
# totpop: https://www2.census.gov/programs-surveys/popest/datasets/2010-2017/national/asrh/nc-est2017-agesex-res.csv

# Code for probabilistic sensitivity analysis (generates distributions of parameters)------------------------------------
if (M > 1) {
        set.seed(314159)
        dirpi17 = matrix(
                data = rnorm(
                        n = M * length(mu_dirpi17),
                        mean = mu_dirpi17,
                        sd = se_dirpi17
                ),
                nrow = length(mu_dirpi17),
                ncol = M
        )
        Npi17j = matrix(data = NA,
                        nrow = 101,
                        ncol = M) #Set up age distribution vector
        for (k in 1:M) {
                Npi17j[1:6, k] =
                        dirpi17[1, k] * totpop[1:6] / sum(totpop[1:6]) #Age 0-5
                Npi17j[7:12, k] =
                        dirpi17[2, k] * totpop[7:12] / sum(totpop[7:12]) #Age 6-11
                Npi17j[13:18, k] =
                        dirpi17[3, k] * totpop[13:18] / sum(totpop[13:18]) #Age 12-17
                Npi17j[19:25, k] =
                        dirpi17[4, k] * totpop[19:25] / sum(totpop[19:25]) #Age 18-24
                Npi17j[26:35, k] =
                        dirpi17[5, k] * totpop[26:35] / sum(totpop[26:35]) #Age 25-34
                Npi17j[36:45, k] =
                        dirpi17[6, k] * totpop[36:45] / sum(totpop[36:45]) #Age 35-44
                Npi17j[46:55, k] =
                        dirpi17[7, k] * totpop[46:55] / sum(totpop[46:55]) #Age 45-54
                Npi17j[56:65, k] =
                        dirpi17[8, k] * totpop[56:65] / sum(totpop[56:65]) #Age 55-64
                Npi17j[66:101, k] = 0
        }
        rho_j = t(t(Npi17j) / apply(Npi17j, 2, sum)) # 101 x M
        sa_rho_j = rho_j #Save for appendix
        theta = lb_theta + (ub_theta - lb_theta) * rbeta(
                n = M,
                shape1 = phi_theta * (mu_theta - lb_theta) / (ub_theta - lb_theta),
                shape2 = phi_theta * (ub_theta - mu_theta) / (ub_theta - lb_theta)
        )
        P =
                rgamma(n = M,
                       shape = shape_P,
                       scale = scale_P) # 1 x M Average $513.16 pmpm for non-APTC consumers
        etaParg_j1 = lb_eta + (ub_eta - lb_eta) * rbeta(
                n = M,
                shape1 = phi_eta1 * (-mu_etaP_j[1] - lb_eta) / (ub_eta - lb_eta),
                shape2 = phi_eta1 * (ub_eta--mu_etaP_j[1]) / (ub_eta - lb_eta)
        )
        etaParg_j2 = lb_eta + (ub_eta - lb_eta) * rbeta(
                n = M,
                shape1 = phi_eta2 * (-mu_etaP_j[2] - lb_eta) / (ub_eta - lb_eta),
                shape2 = phi_eta2 * (ub_eta--mu_etaP_j[2]) / (ub_eta - lb_eta)
        )
        etaParg_j3 = lb_eta + (ub_eta - lb_eta) * rbeta(
                n = M,
                shape1 = phi_eta3 * (-mu_etaP_j[3] - lb_eta) / (ub_eta - lb_eta),
                shape2 = phi_eta3 * (ub_eta--mu_etaP_j[3]) / (ub_eta - lb_eta)
        )
        etaParg_j = rbind(-etaParg_j1, -etaParg_j2, -etaParg_j3) # 3 x M
        etaP_j = rbind(
                matrix(
                        rep((etaParg_j[1,] + etaParg_j[2,]) / 2, 18),
                        nrow = 18,
                        ncol = M,
                        byrow = TRUE
                ),
                matrix(
                        rep(etaParg_j[1,], 17),
                        nrow = 17,
                        ncol = M,
                        byrow = TRUE
                ),
                matrix(
                        rep(etaParg_j[2,], 20),
                        nrow = 20,
                        ncol = M,
                        byrow = TRUE
                ),
                matrix(
                        rep(etaParg_j[3,], 10),
                        nrow = 10,
                        ncol = M,
                        byrow = TRUE
                ),
                matrix(
                        rep(0, 36 * M),
                        nrow = 36,
                        ncol = M,
                        byrow = TRUE
                )
        )
        NNT = rgamma(n = M,
                     shape = shape_NNT,
                     scale = scale_NNT)
        alpha = rbeta(
                n = M,
                shape1 = mu_alpha * phi_alpha,
                shape2 = (1 - mu_alpha) * phi_alpha
        )
        QI_j = t(pmin(boot.predict.out, 1)) #Bootstrapped from NHMS
        sa_QI_j = QI_j #Save for Appendix Table
} else {
        #Code for base case analysis --generates base case parameters---------------------------------------------------------
        dirpi17 = mu_dirpi17
        Npi17j = numeric(length = 101)
        Npi17j[1:6] = dirpi17[1] * totpop[1:6] / sum(totpop[1:6]) #Age 0-5
        Npi17j[7:12] = dirpi17[2] * totpop[7:12] / sum(totpop[7:12]) #Age 6-11
        Npi17j[13:18] = dirpi17[3] * totpop[13:18] / sum(totpop[13:18]) #Age 12-17
        Npi17j[19:25] = dirpi17[4] * totpop[19:25] / sum(totpop[19:25]) #Age 18-24
        Npi17j[26:35] = dirpi17[5] * totpop[26:35] / sum(totpop[26:35]) #Age 25-34
        Npi17j[36:45] = dirpi17[6] * totpop[36:45] / sum(totpop[36:45]) #Age 35-44
        Npi17j[46:55] = dirpi17[7] * totpop[46:55] / sum(totpop[46:55]) #Age 45-54
        Npi17j[56:65] = dirpi17[8] * totpop[56:65] / sum(totpop[56:65]) #Age 55-64
        Npi17j[66:101] = 0
        rho_j = Npi17j / sum(Npi17j) #Age distribution of direct private insurance
        bc_rho_j = rho_j #Save for base case in appendix table
        theta = mu_theta #Proportion of premium passed through to consumers
        P = mu_P #Average annual premium for non-APTC consumers
        etaP_j = c(
                rep((mu_etaP_j[1] + mu_etaP_j[2]) / 2, 18),
                rep(mu_etaP_j[1], 17),
                rep(mu_etaP_j[2], 20),
                rep(mu_etaP_j[3], 10),
                rep(0, 36)
        ) #Elasticity by age from Saltzman 2019
        NNT = mu_NNT #Number needed to lose insurance to result in 1 death
        alpha = mu_alpha #Proportion of disability amenable to health care
        QI_j = mu_QI_j #Age-related HRQoL from NHMS
        bc_QI_j = QI_j #Save for base case in appendix table
}

#Main analysis routine--------------------------------------------------------------------------------------------------
psurv = 1 - phi_j #Survival probability from national life tables
n = 100000 #Number of enrollees / set equal to 1 for probability
deltaC = 10000000 #Cost of new technology / set equal to 1 for marginal effect, equal to 1000000 for interpretability
deltaP = theta * deltaC / n #101 x 1 #Change in premium assuming theta passthrough
r = 0.03 #Discount rate
dv = 1 / (1 + r) ^ seq(1, 101, 1) #Discount vector (1 per year of age into the future)

QAYLL = numeric(length = M)
dQALY = numeric(length = M)
RR = numeric(length = M)
sumNUj = numeric(length = M)
sumNDj = numeric(length = M)
for (i in 1:M) {
        pdp = deltaP[i] / P[i] #%Delta premium assuming theta passthrough of deltaC
        if (M > 1) {
                NUj = -n * rho_j[, i] * etaP_j[, i] * pdp #Number disenrolled per age group
        } else {
                NUj = -n * rho_j * etaP_j * pdp
        }
        RR[i] = max(1 + (sum(NUj[19:65]) / NNT[i]) / sum(phi_j[19:65] * NUj[19:65]), 1) #Rel. risk of mortality (assume >= 1)
        if (M > 1) {
                MI_j = 1 - QI_j[, i] #Morbidity in the insured. Note: QI_j is observed
        } else {
                MI_j = 1 - QI_j
        }
        MU_j = MI_j / ((1 - alpha[i]) + alpha[i] / RR[i])  #Morbidity in the uninsured
        NDj = c((RR[i] - 1) * phi_j[1:18] * NUj[1:18], (RR[i] - 1) * phi_j[19:65] * NUj[19:65], rep(0, 36)) #Number of deaths by age j
        qmat = matrix(0, 101, 101)
        for (j in 1:101) {
                # The sum of each row[j] represents discounted QALYs (insured) lost for each death at age j
                if (M > 1) {
                        qmat[j, j:101] = dv[1:(101 - (j - 1))] * QI_j[j:101, i] * cumprod(psurv[j:101])
                }
                else {
                        qmat[j, j:101] = dv[1:(101 - (j - 1))] * QI_j[j:101] * cumprod(psurv[j:101])
                }
                # Half-cycle correction assuming death uniform through year
                qmat[j, j] = .5 * qmat[j, j]
        }
        QAYLL[i] = sum(apply(X = qmat, MARGIN = 1, sum) * NDj) #Total discounted QALYs lost due to premium-induced disenrollment deaths
        dQALY[i] = sum(NUj * (MU_j - MI_j)) #Disutility increase among newly uninsured for one year
        sumNUj[i] = sum(NUj)
        sumNDj[i] = sum(NDj)
}
Threshold_QALYs = deltaC / (QAYLL + dQALY) #Threshold

if (M == 1) {
        print("Base Case")
        print("Drop Coverage")
        print(sumNUj)
        print("Deaths")
        print(sumNDj)
        print("QALYs lost from mortality")
        print(QAYLL)
        print("QALYs lost from morbidity")
        print(dQALY)
        print("Total QALYs lost")
        print(QAYLL + dQALY)
        bc_Threshold = Threshold_QALYs
        print("Threshold")
        print(bc_Threshold)
        save(
                list = c("bc_Threshold", "bc_rho_j", "bc_QI_j"),
                file = paste("output/threshold_", vstr, "_bc_parms.Rdata", sep = "")
        )
}

#Sensitivity Analysis (only when M > 1)

if (M > 1) {
        load(file = paste("output/threshold_", vstr, "_bc_parms.Rdata", sep = ""))
        print("95% UI:")
        print("Drop Coverage")
        print(quantile(sumNUj, c(.025, .975)))
        print("Deaths")
        print(quantile(sumNDj, c(.025, .975)))
        print("QALYs lost from mortality")
        print(quantile(QAYLL, c(.025, .975)))
        print("QALYs lost from morbidity")
        print(quantile(dQALY, c(.025, .975)))
        print("Total QALYs lost")
        print(quantile(QAYLL + dQALY, c(.025, .975)))
        print("Threshold")
        print(quantile(Threshold_QALYs, c(.025, .975)))
        h = hist(
                Threshold_QALYs,
                breaks = 50,
                main = "",
                xlim = c(0, 350000),
                xlab = "Threshold Value (2019 US$/QALY)",
                ylim = c(0, 7000)
        ) #Figure 1
        cuts = cut(
                h$breaks,
                c(-Inf, 100000, 150000, Inf),
                include.lowest = TRUE,
                right = FALSE
        )
        pdf(
                file = paste("output/PSA_threshold_", vstr, ".pdf", sep = ""),
                width = 8,
                height = 6
        )
        hist(
                Threshold_QALYs,
                breaks = 50,
                main = "",
                xlim = c(0, 350000),
                xlab = "Threshold Value (2019 US$/QALY)",
                col = c("blue", "gray", "orange")[cuts],
                ylim = c(0, 7000)
        )
        segments(
                x0 = 100000,
                y0 = 0,
                x1 = 100000,
                y1 = 6500,
                lty = 1,
                lwd = 2
        )
        arrows(
                x0 = 100000,
                y0 = 3000,
                x1 = 0,
                y1 = 3000,
                code = 2,
                lty = 1
        )
        segments(
                x0 = 150000,
                y0 = 0,
                x1 = 150000,
                y1 = 6500,
                lty = 1,
                lwd = 2
        )
        arrows(
                x0 = 150000,
                y0 = 3000,
                x1 = 325000,
                y1 = 3000,
                code = 2,
                lty = 1
        )
        segments(
                x0 = bc_Threshold,
                y0 = 0,
                x1 = bc_Threshold,
                y1 = 6500,
                lty = 2,
                lwd = 3
        )

        arrows(
                x0 = quantile(Threshold_QALYs, .025),
                y0 = 1000,
                x1 = quantile(Threshold_QALYs, .975),
                y1 = 1000,
                code = 3,
                angle = 90,
                lwd = 3
        )
        text(
                x = bc_Threshold,
                y = 6900,
                labels = paste("Base Case: \n $", round(bc_Threshold, -3), "/QALY", sep = ""),
                cex = 1
        )
        text(
                x = quantile(Threshold_QALYs, .975) + 5000,
                y = 1000,
                labels = paste(
                        "(95% UI: $",
                        round(quantile(Threshold_QALYs, .025),-3),
                        "-$",
                        round(quantile(Threshold_QALYs, .975),-3),
                        "/QALY)",
                        sep = ""
                ),
                adj = c(0, NA),
                cex = 1
        )
        text(
                x = 55000,
                y = 2600,
                adj = c(1, NA),
                labels = paste("Pr = ", round(
                        100 * sum(Threshold_QALYs < 100000) / M, 0
                ), "%", sep = ""),
                cex = 1
        )
        text(
                x = 270000,
                y = 2600,
                adj = c(0, NA),
                labels = paste("Pr = ", round(
                        100 * sum(Threshold_QALYs > 150000) / M, 0
                ), "%", sep = ""),
                cex = 1
        )
        dev.off()

        #PSA Parameters---------------------------------------------------------------------------------------------------------

        nul_ = write.table(
                x = cbind(
                        c("theta",
                          "P",
                          "NNT",
                          "alpha",
                          "eta1",
                          "eta2",
                          "eta3"),
                        rbind(
                                round(quantile(theta, c(.025, .975)), 2),
                                round(quantile(P, c(.025, .975)), 0),
                                round(quantile(NNT, c(.025, .975)), 1),
                                round(quantile(alpha, c(.025, .975)), 3),
                                round(quantile(etaParg_j[1, ], c(.025, .975)), 2),
                                round(quantile(etaParg_j[2, ], c(.025, .975)), 2),
                                round(quantile(etaParg_j[3, ], c(.025, .975)), 2)
                        )
                ),
                file = paste("output/table1_psaparams_", vstr, ".txt", sep = "") #95% UI for parameters in table 1
        )

        nul_ = write.table(
                cbind(
                        round(bc_rho_j, 4),
                        t(as.matrix(round(
                                apply(sa_rho_j, MARGIN = 1, function(x) {
                                        quantile(x, c(.025, .975))
                                }), 4
                        ))),
                        round(phi_j, 4),
                        round(bc_QI_j, 4),
                        t(as.matrix(round(
                                apply(sa_QI_j, MARGIN = 1, function(x) {
                                        quantile(x, c(.025, .975))
                                }), 4
                        )))
                ),
                file = paste(
                        "output/appendix_table1_psaparams_",
                        vstr,
                        ".txt",
                        sep = ""
                ) #95% UI for parameters in appendix table 1
        )

        #Partial Probabilistic SA Parameters------------------------------------------------------------------------------------

        SA_grid = rbind(
                expand.grid(
                        SA_theta = as.numeric(quantile(theta, seq(0, 1, .02))),
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = as.numeric(quantile(P, seq(0, 1, .02))),
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = as.numeric(quantile(NNT, seq(0, 1, .02))),
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = as.numeric(quantile(alpha, seq(0, 1, .02))),
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = as.numeric(quantile(etaParg_j[1, ], seq(0, 1, .02))),
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = as.numeric(quantile(etaParg_j[2, ], seq(0, 1, .02))),
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = as.numeric(quantile(etaParg_j[3, ], seq(0, 1, .02))),
                        SA_r = .03
                ),
                expand.grid(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = etaParg_j[3],
                        SA_r = seq(0, .06, .005)
                )
        )

        #Conduct Partial Probabilistic Sensitivity Analysis-------------------------------------------------------------------------------------
        SA_grid = cbind(SA_grid,
                        SA_etaP_0 = (SA_grid$SA_etaP_1 + SA_grid$SA_etaP_2) / 2)

        dirpi17 = 1000 * c(2358, 2596, 3163, 4327, 5967, 5138, 6675, 7401) #2017 Estimated population covered by direct purchase insurance ages 0-5, 6-11, 12-17, 18-24, 25-34, 35-44, 45-54, 55-64 from Table HI01, 2018 CPS ASES https://www.census.gov/data/tables/time-series/demo/income-poverty/cps-hi/hi-01.html
        Npi17j = numeric(length = 101)
        Npi17j[1:6] = dirpi17[1] * totpop[1:6] / sum(totpop[1:6]) #Age 0-5
        Npi17j[7:12] = dirpi17[2] * totpop[7:12] / sum(totpop[7:12]) #Age 6-11
        Npi17j[13:18] = dirpi17[3] * totpop[13:18] / sum(totpop[13:18]) #Age 12-17
        Npi17j[19:25] = dirpi17[4] * totpop[19:25] / sum(totpop[19:25]) #Age 18-24
        Npi17j[26:35] = dirpi17[5] * totpop[26:35] / sum(totpop[26:35]) #Age 25-34
        Npi17j[36:45] = dirpi17[6] * totpop[36:45] / sum(totpop[36:45]) #Age 35-44
        Npi17j[46:55] = dirpi17[7] * totpop[46:55] / sum(totpop[46:55]) #Age 45-54
        Npi17j[56:65] = dirpi17[8] * totpop[56:65] / sum(totpop[56:65]) #Age 55-64
        Npi17j[66:101] = 0
        rho_j = Npi17j / sum(Npi17j) #Age distribution of direct private insurance
        QI_j = mu_QI_j
        psurv = 1 - phi_j

        SA_Threshold_QALYs = numeric()
        SA_RR = numeric()

        for (k in 1:dim(SA_grid)[1]) {
                etaP_j = c(
                        rep(SA_grid$SA_etaP_0[k], 18),
                        rep(SA_grid$SA_etaP_1[k], 17),
                        rep(SA_grid$SA_etaP_2[k], 20),
                        rep(SA_grid$SA_etaP_3[k], 10),
                        rep(0, 36)
                )
                n = 1
                deltaC = 1
                deltaP = SA_grid$SA_theta[k] * deltaC / n #101 x 1
                dv = 1 / (1 + SA_grid$SA_r[k]) ^ seq(1, 101, 1) #Discount vector (1 per year of age into the future)
                pdp = deltaP / SA_grid$SA_P[k] #%Delta premium assuming theta passthrough of deltaC
                NUj = -n * rho_j * etaP_j * pdp #Number disenrolled per age group
                RR = max(1 + (sum(NUj[19:65]) / SA_grid$SA_NNT[k]) / sum(phi_j[19:65] * NUj[19:65]),
                         1) #Note: floor at 1 - i.e., insurance cannot result in worse mortality
                SA_RR[k] = RR
                MI_j = 1 - QI_j #Morbidity in the insured. Note: QI_j is observed
                MU_j = MI_j / ((1 - SA_grid$SA_alpha[k]) + SA_grid$SA_alpha[k] / RR) #Morbidity in the uninsured
                NDj = c((RR - 1) * phi_j[1:18] * NUj[1:18],
                        (RR - 1) * phi_j[19:65] * NUj[19:65],
                        rep(0, 36)) #Number of deaths by age j
                qmat = matrix(0, 101, 101)
                for (j in 1:101) {
                        # The sum of each row[j] represents discounted QALYs (insured) lost for each death at age j
                        qmat[j, j:101] = dv[1:(101 - (j - 1))] * QI_j[j:101] * cumprod(psurv[j:101])
                        # Credit back assuming death occurs halfway into the year
                        qmat[j, j] = .5 * qmat[j, j]
                }
                QAYLL = sum(apply(X = qmat, MARGIN = 1, sum) * NDj) #Total discounted QALYs lost due to premium-induced disenrollment
                dQALY = sum(NUj * (MU_j - MI_j)) #Disutility increase among newly uninsured for one year
                SA_Threshold_QALYs[k] = deltaC / (QAYLL + dQALY)
        }

        SA_Results = cbind(SA_Threshold_QALYs, SA_grid)

        #Plot Tornado+ diagram--------------------------------------------------------------------------------------------------

        pdf(
                file = paste("output/tornado_", vstr, ".pdf", sep = ""),
                width = 6,
                height = 10
        )
        layout(
                mat = matrix(
                        1:7,
                        nrow = 7,
                        ncol = 1,
                        byrow = TRUE
                ),
                widths = c(1),
                heights = rep(1, 7)
        )
        par(mar = c(.5, 5.1, .5, 1))

        plot(
                x = SA_Results$SA_Threshold_QALYs[103:153],
                y = SA_Results$SA_NNT[103:153],
                ylim = range(SA_Results$SA_NNT),
                xlim = c(0, 300000),
                type = "b",
                ann = TRUE,
                ylab = "NNT",
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )

        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)

        segments(
                x0 = -10000,
                y0 = mu_NNT,
                x1 = bc_Threshold,
                y1 = mu_NNT,
                lty = 2
        )

        par(mar = c(.5, 5.1, 0, 1))

        plot(
                x = SA_Results$SA_Threshold_QALYs[205:255],
                y = SA_Results$SA_etaP_1[205:255],
                ylim = range(SA_Results$SA_etaP_1),
                xlim = c(0, 300000),
                type = "b",
                ann = TRUE,
                ylab = expression(paste(eta ^ P, "18-34")),
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )
        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)
        segments(
                x0 = -10000,
                y0 = mu_etaP_j[1],
                x1 = bc_Threshold,
                y1 = mu_etaP_j[1],
                lty = 2
        )

        plot(
                x = SA_Results$SA_Threshold_QALYs[256:306],
                y = SA_Results$SA_etaP_2[256:306],
                ylim = range(SA_Results$SA_etaP_2),
                xlim = c(0, 300000),
                type = "b",
                ann = TRUE,
                ylab = expression(paste(eta ^ P, "35-54")),
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )
        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)
        segments(
                x0 = -10000,
                y0 = mu_etaP_j[2],
                x1 = bc_Threshold,
                y1 = mu_etaP_j[2],
                lty = 2
        )

        plot(
                x = SA_Results$SA_Threshold_QALYs[1:51],
                y = SA_Results$SA_theta[1:51],
                ylim = range(SA_Results$SA_theta),
                xlim = c(0, 300000),
                type = "b",
                ann = TRUE,
                ylab = expression(theta),
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )
        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)
        segments(
                x0 = -10000,
                y0 = mu_theta,
                x1 = bc_Threshold,
                y1 = mu_theta,
                lty = 2
        )

        plot(
                x = SA_Results$SA_Threshold_QALYs[52:102],
                y = SA_Results$SA_P[52:102],
                ylim = range(SA_Results$SA_P),
                xlim = c(0, 300000),
                type = "b",
                ann = TRUE,
                ylab = "P",
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )
        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)
        segments(
                x0 = -10000,
                y0 = mu_P,
                x1 = bc_Threshold,
                y1 = mu_P,
                lty = 2
        )

        plot(
                x = SA_Results$SA_Threshold_QALYs[154:204],
                y = SA_Results$SA_alpha[154:204],
                ylim = range(SA_Results$SA_alpha),
                xlim = c(0, 300000),
                type = "b",
                ann = TRUE,
                ylab = expression(alpha),
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )
        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)
        segments(
                x0 = -10000,
                y0 = mu_alpha,
                x1 = bc_Threshold,
                y1 = mu_alpha,
                lty = 2
        )

        par(mar = c(2.1, 5.1, 0, 1))

        plot(
                x = SA_Results$SA_Threshold_QALYs[307:357],
                y = SA_Results$SA_etaP_3[307:357],
                ylim = range(SA_Results$SA_etaP_3),
                xlim = c(0, 300000),
                ann = TRUE,
                ylab = expression(paste(eta ^ P, "55-64")),
                cex.axis = 0.8,
                cex.lab = 1.2,
                las = 1
        )
        abline(v = bc_Threshold, lty = 2)
        abline(v = 100000, lty = 1)
        abline(v = 150000, lty = 1)
        segments(
                x0 = -10000,
                y0 = mu_etaP_j[3],
                x1 = bc_Threshold,
                y1 = mu_etaP_j[3],
                lty = 2
        )
        dev.off()

        # Total probabilistic sensitivity analysis --------------------------------------------------------------------------------------------------------------------------------------------------------

        print("Rel to 50000/QALY")
        print(sum((Threshold_QALYs >= (
                100000 - 5000
        )) &
                (Threshold_QALYs <= (
                        100000 + 5000
                ))))
        print(sum((Threshold_QALYs >= (50000 - 5000)) &
                          (Threshold_QALYs <= (50000 + 5000))))
        print(sum((Threshold_QALYs >= (
                100000 - 5000
        )) &
                (Threshold_QALYs <= (
                        100000 + 5000
                ))) /
                sum((Threshold_QALYs >= (50000 - 5000)) &
                            (Threshold_QALYs <= (50000 + 5000))))

        print("Rel to 300000/QALY")
        print(sum((Threshold_QALYs >= (
                100000 - 5000
        )) &
                (Threshold_QALYs <= (
                        100000 + 5000
                ))))
        print(sum((Threshold_QALYs >= (
                300000 - 5000
        )) &
                (Threshold_QALYs <= (
                        300000 + 5000
                ))))
        print(sum((Threshold_QALYs >= (
                100000 - 5000
        )) &
                (Threshold_QALYs <= (
                        100000 + 5000
                ))) /
                sum((Threshold_QALYs >= (
                        300000 - 5000
                )) &
                        (Threshold_QALYs <= (
                                300000 + 5000
                        ))))

        print("P(>150000)")
        print(sum(Threshold_QALYs > 150000))

        print("P(<100000)")
        print(sum(Threshold_QALYs < 100000))


        # Table: Key Input Values and 1-Way Sensitivity Analysis Results ----------------------------------------------------------------------------------------------------------------------------------

        bcpars = list(
                SA_theta = mu_theta,
                SA_P = mu_P,
                SA_NNT = mu_NNT,
                SA_alpha = mu_alpha,
                SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) / 2,
                SA_etaP_1 = mu_etaP_j[1],
                SA_etaP_2 = mu_etaP_j[2],
                SA_etaP_3 = mu_etaP_j[3],
                SA_r = .03
        )

        dirpi17 = 1000 * c(2358, 2596, 3163, 4327, 5967, 5138, 6675, 7401) #2017 Estimated population covered by direct purchase insurance ages 0-5, 6-11, 12-17, 18-24, 25-34, 35-44, 45-54, 55-64 from Table HI01, 2018 CPS ASES https://www.census.gov/data/tables/time-series/demo/income-poverty/cps-hi/hi-01.html
        Npi17j = numeric(length = 101)
        Npi17j[1:6] = dirpi17[1] * totpop[1:6] / sum(totpop[1:6]) #Age 0-5
        Npi17j[7:12] = dirpi17[2] * totpop[7:12] / sum(totpop[7:12]) #Age 6-11
        Npi17j[13:18] = dirpi17[3] * totpop[13:18] / sum(totpop[13:18]) #Age 12-17
        Npi17j[19:25] = dirpi17[4] * totpop[19:25] / sum(totpop[19:25]) #Age 18-24
        Npi17j[26:35] = dirpi17[5] * totpop[26:35] / sum(totpop[26:35]) #Age 25-34
        Npi17j[36:45] = dirpi17[6] * totpop[36:45] / sum(totpop[36:45]) #Age 35-44
        Npi17j[46:55] = dirpi17[7] * totpop[46:55] / sum(totpop[46:55]) #Age 45-54
        Npi17j[56:65] = dirpi17[8] * totpop[56:65] / sum(totpop[56:65]) #Age 55-64
        Npi17j[66:101] = 0
        rho_j = Npi17j / sum(Npi17j) #Age distribution of direct private insurance
        QI_j = mu_QI_j
        psurv = 1 - phi_j

        findthreshold = function(pars) {
                etaP_j = c(
                        rep(pars$SA_etaP_0, 18),
                        rep(pars$SA_etaP_1, 17),
                        rep(pars$SA_etaP_2, 20),
                        rep(pars$SA_etaP_3, 10),
                        rep(0, 36)
                )
                n = 1
                deltaC = 1
                deltaP = pars$SA_theta * deltaC / n #101 x 1
                dv = 1 / (1 + pars$SA_r) ^ seq(1, 101, 1) #Discount vector (1 per year of age into the future)
                pdp = deltaP / pars$SA_P #%Delta premium assuming theta passthrough of deltaC
                NUj = -n * rho_j * etaP_j * pdp #Number disenrolled per age group
                RR = max(1 + (sum(NUj[19:65]) / pars$SA_NNT) / sum(phi_j[19:65] * NUj[19:65]), 1) #Note: floor at 1 - i.e., insurance cannot result in worse mortality
                SA_RR = RR
                MI_j = 1 - QI_j #Morbidity in the insured. Note: QI_j is observed
                MU_j = MI_j / ((1 - pars$SA_alpha) + pars$SA_alpha / RR) #Morbidity in the uninsured
                NDj = c((RR - 1) * phi_j[1:18] * NUj[1:18],
                        (RR - 1) * phi_j[19:65] * NUj[19:65],
                        rep(0, 36)) #Number of deaths by age j
                qmat = matrix(0, 101, 101)
                for (j in 1:101) {
                        # The sum of each row[j] represents discounted QALYs (insured) lost for each death at age j
                        qmat[j, j:101] = dv[1:(101 - (j - 1))] * QI_j[j:101] * cumprod(psurv[j:101])
                        # Credit back assuming death occurs halfway into the year
                        qmat[j, j] = .5 * qmat[j, j]
                }
                QAYLL = sum(apply(X = qmat, MARGIN = 1, sum) * NDj) #Total discounted QALYs lost due to premium-induced disenrollment
                dQALY = sum(NUj * (MU_j - MI_j)) #Disutility increase among newly uninsured for one year
                SA_Threshold_QALYs = deltaC / (QAYLL + dQALY)
                SA_Threshold_QALYs
        }

        # Table Input Values -> (< 100K/QALY or >150K/QALY) ---------------------------------------------------------------------------------------------------------------------------------------------

        #NNT
        NNT = rgamma(n = M,
                     shape = shape_NNT,
                     scale = scale_NNT)

        psa_NNT = function(par_NNT, target) {
                pars_ = bcpars
                pars_$SA_NNT = par_NNT
                (findthreshold(pars_) - target)
        }
        NNT_100000 = uniroot(f = psa_NNT,
                             interval = c(10, 10000),
                             target = 100000)$root
        NNT_150000 = uniroot(f = psa_NNT,
                             interval = c(10, 10000),
                             target = 150000)$root
        print(c(
                NNT_100000,
                pgamma(
                        q = NNT_100000,
                        shape = shape_NNT,
                        scale = scale_NNT
                )
        ))
        print(c(
                NNT_150000,
                1 - pgamma(
                        q = NNT_150000,
                        shape = shape_NNT,
                        scale = scale_NNT
                )
        ))

        #Elasticity 18-34
        psa_etaP_2 = function(par_etaP_2, target) {
                pars_ = bcpars
                pars_$SA_etaP_2 = par_etaP_2
                pars_$SA_etaP_0 = (mu_etaP_j[1] + par_etaP_2) / 2
                (findthreshold(pars_) - target)
        }
        etaP_2_100000 = uniroot(f = psa_etaP_2,
                                interval = c(-3, -.01),
                                target = 100000)$root
        etaP_2_150000 = uniroot(f = psa_etaP_2,
                                interval = c(-3, -.01),
                                target = 150000)$root
        print(c(
                etaP_2_100000,
                1 - pbeta(
                        q = (-etaP_2_100000 - lb_eta) / (ub_eta - lb_eta),
                        shape1 = phi_eta2 * (-mu_etaP_j[2] - lb_eta) / (ub_eta - lb_eta),
                        shape2 = phi_eta2 * (ub_eta--mu_etaP_j[2]) / (ub_eta - lb_eta)
                )
        ))
        print(c(
                etaP_2_150000,
                pbeta(
                        q = (-etaP_2_150000 - lb_eta) / (ub_eta - lb_eta),
                        shape1 = phi_eta2 * (-mu_etaP_j[2] - lb_eta) / (ub_eta - lb_eta),
                        shape2 = phi_eta2 * (ub_eta--mu_etaP_j[2]) / (ub_eta - lb_eta)
                )
        ))

        #Elasticity 35-54
        psa_etaP_1 = function(par_etaP_1, target) {
                pars_ = bcpars
                pars_$SA_etaP_1 = par_etaP_1
                pars_$SA_etaP_0 = (par_etaP_1 + mu_etaP_j[2]) / 2
                (findthreshold(pars_) - target)
        }
        etaP_1_100000 = uniroot(f = psa_etaP_1,
                                interval = c(-3, -.01),
                                target = 100000)$root
        etaP_1_150000 = uniroot(f = psa_etaP_1,
                                interval = c(-3, -.01),
                                target = 150000)$root
        print(c(
                etaP_1_100000,
                1 - pbeta(
                        q = (-etaP_1_100000 - lb_eta) / (ub_eta - lb_eta),
                        shape1 = phi_eta1 * (-mu_etaP_j[1] - lb_eta) / (ub_eta - lb_eta),
                        shape2 = phi_eta1 * (ub_eta--mu_etaP_j[1]) / (ub_eta - lb_eta)
                )
        ))
        print(c(
                etaP_1_150000,
                pbeta(
                        q = (-etaP_1_150000 - lb_eta) / (ub_eta - lb_eta),
                        shape1 = phi_eta1 * (-mu_etaP_j[1] - lb_eta) / (ub_eta - lb_eta),
                        shape2 = phi_eta1 * (ub_eta--mu_etaP_j[1]) / (ub_eta - lb_eta)
                )
        ))

        #Premium passthrough
        psa_theta = function(par_theta, target) {
                pars_ = bcpars
                pars_$SA_theta = par_theta
                (findthreshold(pars_) - target)
        }
        theta_100000 = uniroot(
                f = psa_theta,
                interval = c(lb_theta, ub_theta),
                target = 100000
        )$root
        theta_150000 = uniroot(
                f = psa_theta,
                interval = c(lb_theta, ub_theta),
                target = 150000
        )$root

        print(c(
                theta_100000,
                1 - pbeta(
                        q = (theta_100000 - lb_theta) / (ub_theta - lb_theta),
                        shape1 = phi_theta * (mu_theta - lb_theta) / (ub_theta - lb_theta),
                        shape2 = phi_theta * (ub_theta - mu_theta) / (ub_theta - lb_theta)
                )
        ))

        print(c(
                theta_150000,
                pbeta(
                        q = (theta_150000 - lb_theta) / (ub_theta - lb_theta),
                        shape1 = phi_theta * (mu_theta - lb_theta) / (ub_theta - lb_theta),
                        shape2 = phi_theta * (ub_theta - mu_theta) / (ub_theta - lb_theta)
                )
        ))

        #Baseline premium
        psa_P = function(par_P, target) {
                pars_ = bcpars
                pars_$SA_P = par_P
                (findthreshold(pars_) - target)
        }
        P_100000 = uniroot(f = psa_P,
                           interval = c(10, 10000),
                           target = 100000)$root
        P_150000 = uniroot(f = psa_P,
                           interval = c(10, 10000),
                           target = 150000)$root
        print(c(
                P_100000,
                pgamma(
                        q = P_100000,
                        shape = shape_P,
                        scale = scale_P
                )
        ))
        print(c(
                P_150000,
                1 - pgamma(
                        q = P_150000,
                        shape = shape_P,
                        scale = scale_P
                )
        ))

        #Morbidity amenable to health care
        psa_alpha = function(par_alpha, target) {
                pars_ = bcpars
                pars_$SA_alpha = par_alpha
                (findthreshold(pars_) - target)
        }

        alpha_100000 = uniroot(f = psa_alpha,
                               interval = c(.01, .99),
                               target = 100000)$root
        # alpha_150000 = uniroot(f = psa_alpha, #Note: function will produce error because there is no value for alpha that causes threshold to exceed $150K/Q
        #                        interval = c(.01, .99),
        #                        target = 150000)$root
        alpha_150000 = NA
        print(c(
                alpha_100000,
                1 - pbeta(
                        q = alpha_100000,
                        shape1 = mu_alpha * phi_alpha,
                        shape2 = (1 - mu_alpha) * phi_alpha
                )
        ))
        print(c(
                alpha_150000,
                pbeta(
                        q = alpha_150000,
                        shape1 = mu_alpha * phi_alpha,
                        shape2 = (1 - mu_alpha) * phi_alpha
                )
        ))

        #Elasticity 55-64
        psa_etaP_3 = function(par_etaP_3, target) {
                pars_ = bcpars
                pars_$SA_etaP_3 = par_etaP_3
                (findthreshold(pars_) - target)
        }
        etaP_3_100000 = uniroot(f = psa_etaP_3,
                                interval = c(-3, -.01),
                                target = 100000)$root
        # etaP_3_150000 = uniroot(f = psa_etaP_3, #Note: function will produce error because there is no value for psa_etaP_3 that causes threshold to exceed $150K/Q
        #                         interval = c(-3, -.01),
        #                         target = 150000)$root
        etaP_3_150000 = NA
        print(c(
                etaP_3_100000,
                1 - pbeta(
                        q = (-etaP_3_100000 - lb_eta) / (ub_eta - lb_eta),
                        shape1 = phi_eta3 * (-mu_etaP_j[3] - lb_eta) / (ub_eta - lb_eta),
                        shape2 = phi_eta3 * (ub_eta--mu_etaP_j[3]) / (ub_eta - lb_eta)
                )
        ))
        print(c(
                etaP_3_150000,
                pbeta(
                        q = (-etaP_3_150000 - lb_eta) / (ub_eta - lb_eta),
                        shape1 = phi_eta3 * (-mu_etaP_j[3] - lb_eta) / (ub_eta - lb_eta),
                        shape2 = phi_eta3 * (ub_eta--mu_etaP_j[3]) / (ub_eta - lb_eta)
                )
        ))

        #  Table Threshold 95% UI, 2019 US$/QALY ----------------------------------------------------------------------------------------------------------------------------------------------------------

        #NNT
        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = 155.9,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = 435.1,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )
        # Elasticity 18-34

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = -2.38,
                        SA_etaP_0 = (-2.38 + mu_etaP_j[2]) / 2,
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = -.62,
                        SA_etaP_0 = (-.62 + mu_etaP_j[2]) / 2,
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        # Elasticity 35-54

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = -1.78,
                        SA_etaP_0 = (mu_etaP_j[1]+-1.78) / 2,
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = -0.43,
                        SA_etaP_0 = (mu_etaP_j[1]+-0.43) / 2,
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        # Passthrough \theta

        findthreshold(
                pars = list(
                        SA_theta = .83,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = 1.17,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        # Baseline premium

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = 5147,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = 7369,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        #Morbidity amenable to health care

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = 0.057,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = 0.155,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = mu_etaP_j[3],
                        SA_r = .03
                )
        )

        #Elasticity 55-64

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = -1.23,
                        SA_r = .03
                )
        )

        findthreshold(
                pars = list(
                        SA_theta = mu_theta,
                        SA_P = mu_P,
                        SA_NNT = mu_NNT,
                        SA_alpha = mu_alpha,
                        SA_etaP_0 = (mu_etaP_j[1] + mu_etaP_j[2]) /
                                2,
                        SA_etaP_1 = mu_etaP_j[1],
                        SA_etaP_2 = mu_etaP_j[2],
                        SA_etaP_3 = -0.28,
                        SA_r = .03
                )
        )

}
