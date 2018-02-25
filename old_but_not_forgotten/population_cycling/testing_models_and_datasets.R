
system.time({
  mm.dddd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mm.data,
                                                         preyBirthFxn=preyBirthDensityDependent,
                                                         preyDeathFxn=preyDeathLV,
                                                         predBirthFxn=predBirthDensityDependent1,
                                                         predDeathFxn=predDeathLV,
                                                         n.parameters.in.functions=c(2,1,2,1))
})

system.time({
  mm.dd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mm.data,
                                                       preyBirthFxn=preyBirthDensityDependent,
                                                       preyDeathFxn=preyDeathLV,
                                                       predBirthFxn=predBirthLV,
                                                       predDeathFxn=predDeathLV,
                                                       n.parameters.in.functions=n.bd.params.dd)
})

system.time({
  lh.lv <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=lh.data,
                                                       preyBirthFxn=preyBirthLV,
                                                       preyDeathFxn=preyDeathLV,
                                                       predBirthFxn=predBirthLV,
                                                       predDeathFxn=predDeathLV,
                                                       n.parameters.in.functions=n.bd.params.lk,
                                                       n.initialization.attempts=1000)
})


system.time({
  ra.dd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=ra.data,
                                                       preyBirthFxn=preyBirthDensityDependent,
                                                       preyDeathFxn=preyDeathLV,
                                                       predBirthFxn=predBirthLV,
                                                       predDeathFxn=predDeathLV,
                                                       n.parameters.in.functions=n.bd.params.dd)
})

system.time({
  ra.ddivlev <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=ra.data,
                                                            preyBirthFxn=preyBirthDensityDependent,
                                                            preyDeathFxn=preyDeathIvlev,
                                                            predBirthFxn=predBirthLV,
                                                            predDeathFxn=predDeathLV,
                                                            n.parameters.in.functions=c(2,2,1,1),
                                                            n.initialization.attempts=4000)
})

system.time({
  ra.dddd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=ra.data,
                                                         preyBirthFxn=preyBirthDensityDependent,
                                                         preyDeathFxn=preyDeathLV,
                                                         predBirthFxn=predBirthDensityDependent2,
                                                         predDeathFxn=predDeathLV,
                                                         n.parameters.in.functions=c(2,1,2,1))
})

system.time({
  ra.dddd1 <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=ra.data,
                                                          preyBirthFxn=preyBirthDensityDependent,
                                                          preyDeathFxn=preyDeathLV,
                                                          predBirthFxn=predBirthDensityDependent1,
                                                          predDeathFxn=predDeathLV,
                                                          n.parameters.in.functions=c(2,1,2,1))
})

system.time({
  ra.ddddWatt <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=ra.data,
                                                             preyBirthFxn=preyBirthDensityDependent,
                                                             preyDeathFxn=preyDeathWatt,
                                                             predBirthFxn=predBirthDensityDependent1,
                                                             predDeathFxn=predDeathLV,
                                                             n.parameters.in.functions=c(2,3,2,1),
                                                             n.initialization.attempts=2000)
})


system.time({
  mr.lv <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mr.data,
                                                       preyBirthFxn=preyBirthLV,
                                                       preyDeathFxn=preyDeathLV,
                                                       predBirthFxn=predBirthLV,
                                                       predDeathFxn=predDeathLV,
                                                       n.parameters.in.functions=c(1,1,1,1),
                                                       n.initialization.attempts = 2000)
})

system.time({
  mr.dd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mr.data,
                                                       preyBirthFxn=preyBirthLV,
                                                       preyDeathFxn=preyDeathLV,
                                                       predBirthFxn=predBirthDensityDependent1,
                                                       predDeathFxn=predDeathLV,
                                                       n.parameters.in.functions=c(1,1,2,1),
                                                       n.initialization.attempts = 2000)
})

system.time({
  mr.dddd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mr.data,
                                                         preyBirthFxn=preyBirthDensityDependent,
                                                         preyDeathFxn=preyDeathLV,
                                                         predBirthFxn=predBirthDensityDependent2,
                                                         predDeathFxn=predDeathLV,
                                                         n.parameters.in.functions=c(2,1,2,1),
                                                         n.initialization.attempts = 2000)
})

system.time({
  mr.dddd1 <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mr.data,
                                                          preyBirthFxn=preyBirthDensityDependent,
                                                          preyDeathFxn=preyDeathLV,
                                                          predBirthFxn=predBirthDensityDependent1,
                                                          predDeathFxn=predDeathLV,
                                                          n.parameters.in.functions=c(2,1,2,1),
                                                          n.initialization.attempts = 2000)
})

system.time({
  mr.dd <- fitGeneralizedLotkaVolteraMaximumLikelihood(data=mr.data,
                                                       preyBirthFxn=preyBirthDensityDependent,
                                                       preyDeathFxn=preyDeathLV,
                                                       predBirthFxn=predBirthLV,
                                                       predDeathFxn=predDeathLV,
                                                       n.parameters.in.functions=c(2,1,1,1))
})

visualizeGeneralizedLVModelFit(mm.data,
                               preyBirthFxn=preyBirthDensityDependent,
                               preyDeathFxn=preyDeathLV,
                               predBirthFxn=predBirthLV,
                               predDeathFxn=predDeathLV,
                               mm.dd,
                               relative.y.limit=2)
