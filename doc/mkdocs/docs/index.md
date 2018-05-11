# Welcome to Π4U

Π4U is a high performance framework for optimization, sampling and Bayesian uncertainty quantification of large scale computational models.

The framework is based on the [TORC](references.md#anchor_torc) task-parallel library for clusters, which is designed to provide unified programming and runtime support for computing platforms that range from single-core systems to hybrid multicore-GPU clusters and heterogenous Grid based supercomputers.




<br><br>


# What Π4U can do for you

1. **Optimize**: given a cost function $F(\vartheta)$ find
	$$
	\vartheta^\star = \mathop{\arg\min}\limits_{\vartheta} F(\vartheta) \,.
	$$

2. **Sample**: given the density of a probability distribution $p_{\vartheta}$ draw samples,
	$$
		\vartheta^{(k)} \sim p_\vartheta, \quad k=1,\ldots,N_s \, .
	$$

3. **Uncertainty Quantification**: given a set of data $d$, the output of the model $f(x;\vartheta)$ a likelihood function $p(d|\vartheta)$ and a prior probablity density $p(\vartheta)$ sample the posterior distribution,
	$$
	p(\vartheta | d) = \frac{p(d | \vartheta) p(\vartheta)}{p(d)}\, .
	$$
The model output $f$ depends on a set of input parameters $x$.

After [installing](installation.md) the software have a look at the [examples](./examples/sampling.md) and learn how you can run your cases.







<br><br><br><br><br><br>

!!! warning
    The software and the documentation page are under continuous development. New pages and new feature will be constantly added.

<!--
# Additional documentation
* Tutorial: [pdf](http://www.cse-lab.ethz.ch/images/software/Pi4Ututorial.pdf)
* Poster about Pi4U: [pdf](http://www.cse-lab.ethz.ch/images/stories/articles/Pi4U-Poster.pdf)
* Presentation at the Europar 2015 conference: [pdf](http://www.cse-lab.ethz.ch/images/stories/Publications/2015/Pi4U.Europar2015.key.pdf)
-->
