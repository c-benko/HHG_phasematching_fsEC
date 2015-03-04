# abs_limit.py
from pylab import *
close('all')
matplotlib.rcParams.update({'font.size': 22})

def abs_limit(Lcoh, Labs, Lmed):
    fac1 = 4 * Labs / (1 + 4 * pi ** 2 * (Labs/Lcoh) ** 2)
    fac2 = 1 + exp(-Lmed / Labs) - 2 * cos(pi * Lmed / Lcoh) * exp(-Lmed / 2 / Labs)
    return  fac1 * fac2 / 4


lmed_span = linspace(0, 15, 1000)

figure(figsize = (10,8))
plot(lmed_span, abs_limit(1,1,lmed_span), label = 'Lcoh = Labs', linewidth = 2)
plot(lmed_span, abs_limit(5,1,lmed_span), label = 'Lcoh = 5 Labs', linewidth = 2)
plot(lmed_span, abs_limit(10,1,lmed_span), label = 'Lcoh = 10 Labs', linewidth = 2)
plot(lmed_span, abs_limit(1000,1,lmed_span), label = 'Lcoh >> Labs', linewidth = 2)
plot(lmed_span, lmed_span ** 2 / 4, label = 'Labs = 0', linewidth = 2)
ylim(0,1)
xlim(0,15)
legend(loc = 4, prop={'size':18})
xlabel('Medium Length')
ylabel('Harmonic Yield')
grid()
show()