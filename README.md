# VaR-and-Greeks-of-Dispersion-Trade-Portfolio
Consider a dispersion trade that involves writing DJX (Dow Jones Index) straddles, and buying straddles on the index components. You will hold the options on the index components in proportion to their underlying stocks’ weights in the index. In a vega-neutral dispersion trade, the sum of the vegas of the index and component options is zero. The inputs of this project include the prices of the 62 options (1 index put, 1 index call, 30 component puts, and 30 component calls), relevant historical data on the returns of the underlying stocks. The outputs include the value-at-risk of the portfolio of options and the various Greek letter risks of the portfolio.