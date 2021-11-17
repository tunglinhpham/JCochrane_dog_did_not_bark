# John Cochrane's "The Dog That Did Not Bark: A Defense of Return Predictability" replicating project
This is a project replicating the result of John Cochrane's famous paper about return's predictability (https://www.jstor.org/stable/40056861).
In this paper, Cochrane found that "the absence of dividend growth predictability gives stronger evidence than does the presence of return predictability".
Also, another finding of Cochrane is that "Long-horizon return forecasts give the same strong evidence".

The dataset used in this project is the annual Center for Research in Security (CRSP) data from 1926 to 2017.

The first table demonstrate the results of the regressions:
The first 3 rows are actual return against past dividend/price ratio, actual risk free return against past dividend/price ratio and dividend growth against past dividend/price ratio.
The last 2 rows are log of return against log of past dividend/price ratio and log of dividend growth against log of past dividend/price ratio.
![Table 1 - Forecasting Regressions](https://github.com/tunglinhpham/JCochrane_dog_did_not_bark/blob/main/Tables%20%26%20Figures/Table%201.png)
The result is quite similar to the paper’s original result.
The coefficients of dividend growth and log dividend growth are negative as expected, as high dividend/price ratio implies low price, which indicates that the dividend in the future would be low.
