# Overview
This repository covers codes and data implemented in the research paper titled ''Do We Need Higher-Order Comoments to Enhance Mean-Variance Portfolios? Evidence from a Simplified Jump Process'' The paper is available on [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3523379).

## Data
The data was collected from two sources. The first one is CRSP, which is based on scraping the S&P 500 constituents. The CRSP data is available [here](https://github.com/simaan84/higher_order_comoment/blob/main/snp_crsp.csv). The second is more common and relies on the book to market portfolios by Fama-French. The second data is available [here](https://raw.githubusercontent.com/simaan84/higher_order_comoment/main/FF_BTM.csv).

## Codes
The repository contains two main codes. The first code illustrates how to construct the final time series for estimation and portfolio formation. This file is named "create_data.R". The second file contains the main code of the research. It includes estimation functions, portfolio optimization, backtesting, and bootstrapping. The main file reproduces all results of the manuscript. 
