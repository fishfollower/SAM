  ----------------------------- -----------------
  catchplot {stockassessment}   R Documentation
  ----------------------------- -----------------

SAM catch plot
--------------

### Description

SAM catch plot

### Usage

    catchplot(fit, obs.show = TRUE, drop = NULL, ...)

### Arguments

+--------------------------------------+--------------------------------------+
| `fit`                                | the object returned from sam.fit     |
+--------------------------------------+--------------------------------------+
| `obs.show`                           | if observations are to be shown also |
+--------------------------------------+--------------------------------------+
| `drop`                               | number of years to be left unplotted |
|                                      | at the end. Default (NULL) is to not |
|                                      | show years at the end with no catch  |
|                                      | information                          |
+--------------------------------------+--------------------------------------+
| `...`                                | extra arguments transferred to plot  |
|                                      | including the following:\            |
|                                      |  `add` logical, plotting is to be    |
|                                      | added on existing plot\              |
|                                      |  `ci` logical, confidence intervals  |
|                                      | should be plotted\                   |
|                                      |  `cicol` color to plot the           |
|                                      | confidence polygon                   |
+--------------------------------------+--------------------------------------+

### Details

Plot of estimated (and optionally observed) total catch in weight

------------------------------------------------------------------------

<div style="text-align: center;">

\[Package *stockassessment* version 0.5.2 [Index](00Index.html)\]

</div>
