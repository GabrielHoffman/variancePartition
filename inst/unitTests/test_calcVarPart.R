

test_calcVarPart = function(){

	library(variancePartition)
	library(lme4)
	library(MASS)

	fit = lm(Reaction ~ Days, sleepstudy)
	calcVarPart(fit)

	fit = glm(Reaction ~ Days, data=sleepstudy)
	calcVarPart(fit)

	fit = lmer(Reaction ~ Days + (1|Subject), data=sleepstudy)
	calcVarPart(fit)

	## Dobson (1990) Page 93: Randomized Controlled Trial :
	counts <- c(18,17,15,20,10,20,25,13,12)
	outcome <- gl(3,1,9)
	treatment <- gl(3,3)
	df = data.frame(treatment, outcome, counts) # showing data
	glm.D93 <- glm(counts ~ outcome + treatment, data = df, family = poisson())
	calcVarPart(glm.D93)

	fitnb <- glm.nb(counts ~ outcome + treatment, data = df)
	calcVarPart(fitnb)

	# Fixed effects models
	######################

	# binomial
	gm1 <- glm(cbind(incidence, size - incidence) ~ period + herd,
	                   data = cbpp, family = binomial)
	calcVarPart(gm1)

	# Poisson
	gm1 <- glm(incidence ~ offset(log(size)) + period + herd,
	                   data = cbpp, family = poisson)
	calcVarPart(gm1)


	# negative binomial mixed model
	gm1 <- glm.nb(incidence ~ offset(log(size)) + period + herd,
	                   data = cbpp)
	calcVarPart(gm1)

	# mixed models
	##############

	# binomial
	gm1 <- glmer(cbind(incidence, size - incidence) ~ (1|period) + (1 | herd),
	                   data = cbpp, family = binomial)
	calcVarPart(gm1)

	# Poisson
	gm1 <- glmer(incidence ~ offset(log(size)) + (1|period) + (1 | herd),
	                   data = cbpp, family = poisson)
	calcVarPart(gm1)


	# negative binomial
	gm1 <- glmer.nb(incidence ~ offset(log(size)) + (1|period) + (1 | herd),
	                   data = cbpp)
	calcVarPart(gm1)

	y = with(cbpp, log(incidence+.5) - log(size))
	gm1 <- lmer(y ~ (1|period) + (1 | herd), data = cbpp)
	calcVarPart(gm1)

	TRUE
}
