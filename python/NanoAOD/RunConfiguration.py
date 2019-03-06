
Runyear = 2017

if Runyear == 2016:
    import Samplelist as Slist
elif Runyear == 2017:
    import Samplelist2017 as Slist

IntLumiyear = {
    2016 : 35.92,
    2017 : 41.53,
    2018 : 59.97
}

IntLumi = IntLumiyear[Runyear]

