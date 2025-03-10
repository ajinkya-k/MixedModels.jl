
# Benchmark Report for _/home/bates/.julia/packages/MixedModels/dn0WY/src/MixedModels.jl_ {#Benchmark-Report-for-/home/bates/.julia/packages/MixedModels/dn0WY/src/MixedModels.jl}

## Job Properties {#Job-Properties}
- Time of benchmark: 2 Oct 2018 - 13:42
  
- Package commit: non gi
  
- Julia commit: 5d4eac
  
- Julia command flags: None
  
- Environment variables: None
  

## Results

Below is a table of this job&#39;s results, obtained by running the benchmarks. The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to index into the BaseBenchmarks suite to retrieve the corresponding benchmarks. The percentages accompanying time and memory values in the below table are noise tolerances. The &quot;true&quot; time/memory value for a given benchmark is expected to fall within this percentage of the reported value. An empty cell means that the value was zero.

|                                                                         ID |                  time |         GC time |          memory |     allocations |
| --------------------------------------------------------------------------:| ---------------------:| ---------------:| ---------------:| ---------------:|
|                              `[&quot;crossed&quot;, &quot;Assay:1+A+B*C+(1 |                 G)+(1 |      H)&quot;]` |   2.943 ms (5%) |                 |
|                           `[&quot;crossed&quot;, &quot;Demand:1+U+V+W+X+(1 |                 G)+(1 |      H)&quot;]` |   2.775 ms (5%) |                 |
|                             `[&quot;crossed&quot;, &quot;InstEval:1+A*I+(1 |                 G)+(1 |      H)&quot;]` |    1.247 s (5%) |      114.131 ms |
|                               `[&quot;crossed&quot;, &quot;InstEval:1+A+(1 |                 G)+(1 |           H)+(1 |      I)&quot;]` |    1.999 s (5%) |
|                               `[&quot;crossed&quot;, &quot;Penicillin:1+(1 |                 G)+(1 |      H)&quot;]` |   2.697 ms (5%) |                 |
|                           `[&quot;crossed&quot;, &quot;ScotsSec:1+A+U+V+(1 |                 G)+(1 |      H)&quot;]` |   4.833 ms (5%) |                 |
|                    `[&quot;crossed&quot;, &quot;dialectNL:1+A+T+U+V+W+X+(1 |                 G)+(1 |           H)+(1 |      I)&quot;]` | 416.892 ms (5%) |
|                           `[&quot;crossed&quot;, &quot;egsingle:1+A+U+V+(1 |                 G)+(1 |      H)&quot;]` |  31.421 ms (5%) |        3.427 ms |
|                                     `[&quot;crossed&quot;, &quot;ml1m:1+(1 |                 G)+(1 |      H)&quot;]` |   36.714 s (5%) |      225.872 ms |
|                            `[&quot;crossed&quot;, &quot;paulsim:1+S+T+U+(1 |                 H)+(1 |      G)&quot;]` |  14.097 ms (5%) |                 |
|                 `[&quot;crossedvector&quot;, &quot;bs10:1+U+V+W+((1+U+V+W) |         G)+((1+U+V+W) |      H)&quot;]` | 165.171 ms (5%) |        3.149 ms |
|                           `[&quot;crossedvector&quot;, &quot;d3:1+U+((1+U) |             G)+((1+U) |       H)+((1+U) |      I)&quot;]` |   49.023 s (5%) |
|                               `[&quot;crossedvector&quot;, &quot;d3:1+U+(1 |                 G)+(1 |           H)+(1 |      I)&quot;]` | 299.348 ms (5%) |
|         `[&quot;crossedvector&quot;, &quot;gb12:1+S+T+U+V+W+X+Z+((1+S+U+W) |         G)+((1+S+T+V) |      H)&quot;]` | 134.101 ms (5%) |                 |
| `[&quot;crossedvector&quot;, &quot;kb07:1+S+T+U+V+W+X+Z+((1+S+T+U+V+W+X+Z) | G)+((1+S+T+U+V+W+X+Z) |      H)&quot;]` |    3.488 s (5%) |       16.508 ms |
|                 `[&quot;crossedvector&quot;, &quot;kb07:1+S+T+U+V+W+X+Z+(1 |             G)+((0+S) |       G)+((0+T) |       G)+((0+U) |       G)+((0+V) |
|                                    `[&quot;nested&quot;, &quot;Animal:1+(1 |                 G)+(1 |      H)&quot;]` |   1.261 ms (5%) |                 |
|                                    `[&quot;nested&quot;, &quot;Chem97:1+(1 |                 G)+(1 |      H)&quot;]` |  58.460 ms (5%) |        6.975 ms |
|                                  `[&quot;nested&quot;, &quot;Chem97:1+U+(1 |                 G)+(1 |      H)&quot;]` |  59.353 ms (5%) |        7.019 ms |
|                                `[&quot;nested&quot;, &quot;Genetics:1+A+(1 |                 G)+(1 |      H)&quot;]` |   2.062 ms (5%) |                 |
|                                    `[&quot;nested&quot;, &quot;Pastes:1+(1 |                 G)+(1 |      H)&quot;]` |   2.298 ms (5%) |                 |
|                                   `[&quot;nested&quot;, &quot;Semi2:1+A+(1 |                 G)+(1 |      H)&quot;]` |   2.309 ms (5%) |                 |
|                         `[&quot;simplescalar&quot;, &quot;Alfalfa:1+A*B+(1 |            G)&quot;]` |   1.210 ms (5%) |                 | 208.80 KiB (1%) |
|                         `[&quot;simplescalar&quot;, &quot;Alfalfa:1+A+B+(1 |            G)&quot;]` |   1.021 ms (5%) |                 | 168.47 KiB (1%) |
|                    `[&quot;simplescalar&quot;, &quot;AvgDailyGain:1+A*U+(1 |            G)&quot;]` |   1.287 ms (5%) |                 | 193.33 KiB (1%) |
|                    `[&quot;simplescalar&quot;, &quot;AvgDailyGain:1+A+U+(1 |            G)&quot;]` |   1.144 ms (5%) |                 | 169.59 KiB (1%) |
|                             `[&quot;simplescalar&quot;, &quot;BIB:1+A*U+(1 |            G)&quot;]` |   1.574 ms (5%) |                 | 222.20 KiB (1%) |
|                             `[&quot;simplescalar&quot;, &quot;BIB:1+A+U+(1 |            G)&quot;]` |   1.171 ms (5%) |                 | 171.31 KiB (1%) |
|                              `[&quot;simplescalar&quot;, &quot;Bond:1+A+(1 |            G)&quot;]` | 958.770 μs (5%) |                 | 141.25 KiB (1%) |
|                     `[&quot;simplescalar&quot;, &quot;Cultivation:1+A*B+(1 |            G)&quot;]` |   1.089 ms (5%) |                 | 173.38 KiB (1%) |
|                       `[&quot;simplescalar&quot;, &quot;Cultivation:1+A+(1 |            G)&quot;]` |   1.138 ms (5%) |                 | 162.14 KiB (1%) |
|                     `[&quot;simplescalar&quot;, &quot;Cultivation:1+A+B+(1 |            G)&quot;]` |   1.147 ms (5%) |                 | 173.47 KiB (1%) |
|                           `[&quot;simplescalar&quot;, &quot;Dyestuff2:1+(1 |            G)&quot;]` | 830.840 μs (5%) |                 | 105.20 KiB (1%) |
|                            `[&quot;simplescalar&quot;, &quot;Dyestuff:1+(1 |            G)&quot;]` | 974.091 μs (5%) |                 | 120.86 KiB (1%) |
|                          `[&quot;simplescalar&quot;, &quot;Exam:1+A*U+B+(1 |            G)&quot;]` |   2.250 ms (5%) |                 |   1.17 MiB (1%) |
|                          `[&quot;simplescalar&quot;, &quot;Exam:1+A+B+U+(1 |            G)&quot;]` |   2.133 ms (5%) |                 |   1.03 MiB (1%) |
|                          `[&quot;simplescalar&quot;, &quot;Gasoline:1+U+(1 |            G)&quot;]` |   1.164 ms (5%) |                 | 162.03 KiB (1%) |
|                       `[&quot;simplescalar&quot;, &quot;Hsb82:1+A+B+C+U+(1 |            G)&quot;]` |   3.048 ms (5%) |                 |   2.12 MiB (1%) |
|                    `[&quot;simplescalar&quot;, &quot;IncBlk:1+A+U+V+W+Z+(1 |            G)&quot;]` |   1.226 ms (5%) |                 | 208.83 KiB (1%) |
|                       `[&quot;simplescalar&quot;, &quot;Mississippi:1+A+(1 |            G)&quot;]` | 980.968 μs (5%) |                 | 145.75 KiB (1%) |
|                              `[&quot;simplescalar&quot;, &quot;PBIB:1+A+(1 |            G)&quot;]` |   1.509 ms (5%) |                 | 234.47 KiB (1%) |
|                                `[&quot;simplescalar&quot;, &quot;Rail:1+(1 |            G)&quot;]` |   1.251 ms (5%) |                 | 151.34 KiB (1%) |
|                   `[&quot;simplescalar&quot;, &quot;Semiconductor:1+A*B+(1 |            G)&quot;]` |   1.313 ms (5%) |                 | 222.95 KiB (1%) |
|            `[&quot;simplescalar&quot;, &quot;TeachingII:1+A+T+U+V+W+X+Z+(1 |            G)&quot;]` |   1.483 ms (5%) |                 | 284.53 KiB (1%) |
|                            `[&quot;simplescalar&quot;, &quot;cake:1+A*B+(1 |            G)&quot;]` |   1.606 ms (5%) |                 | 412.83 KiB (1%) |
|                         `[&quot;simplescalar&quot;, &quot;ergoStool:1+A+(1 |            G)&quot;]` |   1.057 ms (5%) |                 | 155.59 KiB (1%) |
|                 `[&quot;singlevector&quot;, &quot;Early:1+U+U&amp;A+((1+U) |            G)&quot;]` |  20.373 ms (5%) |                 |   3.47 MiB (1%) |
|                        `[&quot;singlevector&quot;, &quot;HR:1+A*U+V+((1+U) |            G)&quot;]` |   5.183 ms (5%) |                 | 915.00 KiB (1%) |
|                        `[&quot;singlevector&quot;, &quot;Oxboys:1+U+((1+U) |            G)&quot;]` |  13.207 ms (5%) |                 |   1.93 MiB (1%) |
|                          `[&quot;singlevector&quot;, &quot;SIMS:1+U+((1+U) |            G)&quot;]` |  61.675 ms (5%) |                 |  12.86 MiB (1%) |
|                        `[&quot;singlevector&quot;, &quot;WWheat:1+U+((1+U) |            G)&quot;]` |   7.311 ms (5%) |                 | 902.31 KiB (1%) |
|                     `[&quot;singlevector&quot;, &quot;Weights:1+A*U+((1+U) |            G)&quot;]` |  18.303 ms (5%) |                 |   3.20 MiB (1%) |
|                    `[&quot;singlevector&quot;, &quot;sleepstudy:1+U+((1+U) |            G)&quot;]` |   4.829 ms (5%) |                 | 797.48 KiB (1%) |
|                        `[&quot;singlevector&quot;, &quot;sleepstudy:1+U+(1 |             G)+((0+U) |      G)&quot;]` |   3.219 ms (5%) |                 |


## Benchmark Group List {#Benchmark-Group-List}

Here&#39;s a list of all the benchmark groups executed by this job:
- `["crossed"]`
  
- `["crossedvector"]`
  
- `["nested"]`
  
- `["simplescalar"]`
  
- `["singlevector"]`
  

## Julia versioninfo {#Julia-versioninfo}

```
Julia Version 1.0.0
Commit 5d4eaca0c9 (2018-08-08 20:58 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
      Ubuntu 18.04.1 LTS
  uname: Linux 4.15.0-36-generic #39-Ubuntu SMP Mon Sep 24 16:19:09 UTC 2018 x86_64 x86_64
  CPU: Intel(R) Core(TM) i5-3570 CPU @ 3.40GHz: 
              speed         user         nice          sys         idle          irq
       #1  1690 MHz     140498 s        134 s      18382 s    1495130 s          0 s
       #2  2513 MHz     131505 s         16 s      18277 s    1504212 s          0 s
       #3  1900 MHz     145131 s        581 s      18892 s    1485409 s          0 s
       #4  1682 MHz     190751 s         38 s      17941 s    1445446 s          0 s
       
  Memory: 15.554645538330078 GB (10502.1171875 MB free)
  Uptime: 16578.0 sec
  Load Avg:  1.4091796875  2.07080078125  1.63037109375
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.0 (ORCJIT, ivybridge)
```

