# John Cochrane's "The Dog That Did Not Bark: A Defense of Return Predictability" replicating project
This is a project replicating the result of John Cochrane's famous paper about return's predictability (https://www.jstor.org/stable/40056861).
In this paper, Cochrane found that "the absence of dividend growth predictability gives stronger evidence than does the presence of return predictability".
Also, another finding of Cochrane is that "Long-horizon return forecasts give the same strong evidence".

The dataset used in this project is the annual Center for Research in Security (CRSP) data from 1926 to 2017.

## Table 1

The first table demonstrates the results of the regressions:
The first 3 rows are actual return (R<sub>t+1</sub>) against past dividend/price ratio (D<sub>t</sub>/P<sub>t</sub>), actual risk free return (R<sub>t+1</sub> - R<sub>t</sub><sup>f</sup>) against past dividend/price ratio (D<sub>t</sub>/P<sub>t</sub>), and dividend growth (D<sub>t+1</sub>/P<sub>t+1</sub>) against past dividend/price ratio (D<sub>t</sub>/P<sub>t</sub>).
The last 2 rows are log of return (r<sub>t+1</sub>) against log of past dividend/price ratio (d<sub>t</sub>/p<sub>t</sub>), and log of dividend growth (Δd<sub>t+1</sub>) against log of past dividend/price ratio (d<sub>t</sub>/p<sub>t</sub>). ***From now on, the smaller letters represent the log of its corresponding capital letter***.
![Table 1 - Forecasting Regressions](https://github.com/tunglinhpham/JCochrane_dog_did_not_bark/blob/main/Tables%20%26%20Figures/Table%201.png)
The result is quite similar to the paper’s original result.
The coefficients of dividend growth and log dividend growth are negative as expected, as high dividend/price ratio implies low price, which indicates that the dividend in the future would be low. Interestingly, the coefficients in the equation of Dividend growth (equation 3) and log of Dividend growth (equation 5) are not statistically significant, which means dividend growth is not forecastable.
