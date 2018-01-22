# Hypernets R Package
Estimate and work with Generalised Hypergeometric Network Ensemnbles

## Installation
To install from a private repo, generate a token from https://github.com/settings/tokens. 
You only need the repo scope.

In R, run the following commands to install the package:
```
install.packages('devtools')
library(devtools)
install_github("sg-dev/hypernets", auth_token = 'your_token_goes_here')

library(hypernets)
```
Best practice is to save your token in env var called GITHUB_PAT.
