# HQCD-P
A package for the Pomeron in Holographic QCD

# Redis dependence
This packages requires [redis](https://redis.io/) running on your computer or in a local accessible one. Redis is an *in memory* database which is used to cache long computations. You can provide the IP of the redis instance in the `init()` function or disable it. However it is highly recommended that you use the cache system provided.
