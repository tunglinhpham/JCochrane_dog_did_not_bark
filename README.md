# John Cochrane's "The Dog That Did Not Bark: A Defense of Return Predictability" replicating project
This is a project replicating the result of John Cochrane's famous paper about return's predictability (https://www.jstor.org/stable/40056861).
In this paper, Cochrane found that "the absence of dividend growth predictability gives stronger evidence than does the presence of return predictability".
Also, another finding of Cochrane is that "Long-horizon return forecasts give the same strong evidence".

The dataset used in this project is the annual Center for Research in Security (CRSP) data from 1926 to 2017.

The first question is: "Are stock returns predictable?". An attempt to answer the question is demonstrated in Table 1:

## Table 1

![Table 1 - Forecasting Regressions](https://github.com/tunglinhpham/JCochrane_dog_did_not_bark/blob/main/Tables%20%26%20Figures/Table%201.png)
R<sub>t+1</sub> is the real return (adjusted by inflation), R<sup>f</sup><sub>t+1</sub> is the real risk-free return, D<sub>t+1</sub>/D<sub>t</sub> is the dividend growth, D<sub>t</sub>/P<sub>t</sub> is past dividend/price ratio. Smaller letters represent the log of its corresponding capital letter.

In term of economic, the Table 1 shows that return **is predictable**: if dividend yield increases by 1%, the price of the stock would increase by about 2% so the return would increase by about 3%. This result contradicts the unpredictability of return: If return is unpredictable, when dividend yield increases by 1%, the price would decrease by 1% so that return be rendered unpredictable (R<sub>t+1</sub> = (D<sub>t+1</sub> + P<sub>t+1</sub>)/P<sub>t</sub>). Also, the sign of the coefficients in the dividend growth regressions is negative as expected: a high high dividend/price ratio implies low price, which indicates that the dividend in the future would be low.

However, in term of statistical significance, the predictability of return is not really strong, with the t-stat of the coefficients of return (or excess return) is slightly above 2. Moreover, several authors pointed out that the regressions about return in Table 1 is biased upward, while the t-stat is biased toward rejection. Therefore, the reliability of those regressions are questionable, or in other word, "return forecastability is dead". Thankfully, Professor John Cochrane provided some stronger tests to see whether the return is forecastable or not (the tests are not about how to better forecast return).

Professor Cochrane's argument is:
>If both returns and dividend growth are unforecastable, then present value logic implies that the price/dividend ratio is constant, which it
obviously is not.
Therefore, the question is not just "Are returns forecastable?" or "Is dividend growth forecastable?", but "*Which* of dividend growth or return is forecastable?", or more specifically, "How much of each?".


