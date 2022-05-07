library(pracma)

##fishing


##move

pop = init()
pop2 = init()
SST.patches = init_SST(years = 84)
pop = spawn(pop)
pop2 = spawn(pop2)
pop = recruit(pop) #s1 ~ 0.002443293
pop2 = recruit2(pop2)
pop
pop2
pop - pop2
mean(pop - pop2)
