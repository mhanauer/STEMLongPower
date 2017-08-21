Show in New WindowClear OutputExpand/Collapse Output
          [,1]      [,2]      [,3]      [,4]      [,5]
[1,] 0.8366600 0.6576610 0.5817954 0.6766035 0.6779079
[2,] 0.6576610 0.8366600 0.6784838 0.6102870 0.5913229
[3,] 0.5817954 0.6784838 0.8366600 0.6145267 0.5205849
[4,] 0.6766035 0.6102870 0.6145267 0.8366600 0.6909007
[5,] 0.6880935 0.5913229 0.5205849 0.6909007 0.8366600
$values
[1] 3.3605311 0.3810294 0.1968399 0.1353238 0.1095760

$vectors
           [,1]       [,2]         [,3]        [,4]       [,5]
[1,] -0.4574110 -0.2250287 -0.433347921  0.69225648  0.3038288
[2,] -0.4488511  0.4032239 -0.540206048 -0.37933281 -0.5070576
[3,] -0.4288251  0.6498819  0.389976473  0.07736302  0.4797034
[4,] -0.4569527 -0.2475851  0.606837949  0.17043823 -0.5649702
[5,] -0.4434086 -0.5505704 -0.007284169 -0.58468227  0.3182357

[1] 1000    4
R Console
 
 
TIME
<dbl>
ITP
<dbl>
ID
<int>
WHITE
<dbl>
1.1	1	4.586286	1	1
2.1	1	3.609644	2	1
3.1	1	2.082326	3	1
4.1	1	3.216616	4	1
5.1	1	3.679696	5	1
6.1	1	1.625839	6	1
6 rows
data.frame
6 x 4
 
 
TIME
<dbl>
ITP
<dbl>
ID
<int>
WHITE
<dbl>
1.1	1	4.586286	1	1
2.1	1	3.609644	2	1
3.1	1	2.082326	3	1
4.1	1	3.216616	4	1
5.1	1	3.679696	5	1
6.1	1	1.625839	6	1
6 rows
Show in New WindowClear OutputExpand/Collapse Output
Linear mixed model fit by REML ['lmerMod']
Formula: ITP ~ TIME + (TIME | ID)
   Data: dataTest

REML criterion at convergence: 1789.1

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-2.5192 -0.6127 -0.0258  0.6064  2.9142 

Random effects:
 Groups   Name        Variance  Std.Dev. Corr 
 ID       (Intercept) 5.862e-01 0.765622      
          TIME        6.729e-05 0.008203 -1.00
 Residual             2.033e-01 0.450880      
Number of obs: 1000, groups:  ID, 200

Fixed effects:
            Estimate Std. Error t value
(Intercept) 3.501103   0.063632   55.02
TIME        0.002476   0.010099    0.25

Correlation of Fixed Effects:
     (Intr)
TIME -0.523
Modify Chunk OptionsRun All Chunks AboveRun Current Chunk
Show in New WindowClear OutputExpand/Collapse Output
Power for predictor 'TIME', (95% confidence interval):
      100.0% (69.15, 100.0)

Test: Kenward Roger (package pbkrtest)
      Effect size for TIME is 0.40

Based on 10 simulations, (0 warnings, 0 errors)
alpha = 0.05, nrow = 1000

Time elapsed: 0 h 0 m 24 s
Modify Chunk OptionsRun All Chunks AboveRun Current Chunk
Show in New WindowClear OutputExpand/Collapse Output
 Show Traceback
Error in getDefaultXname(fit) : Couldn't automatically determine a default fixed effect for this model.
