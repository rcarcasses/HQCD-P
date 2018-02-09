init()
# compute this with 400 cheb points
p <- HQCDP()
p <- addKernel(p, potential = UJgTest,
               numReg = 4,
               comment = 'Leading twist gluon sector',
               kernelName = 'gluon', # this has to be unique: is used to name the couplings and the kernel
               optimPars = c(invls = 1/0.153, a = -4.35, b = 1.41, c = 0.626, d = -0.117))
p <- addProcessObservable(p, F2())
p <- addProcessObservable(p, DVCSDSigma())
fit(p)
