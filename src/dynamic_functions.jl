const g = 5000
const SigmaRes = 6*mPi^2
const ϵ = 5e-5

```
	Two-particle interaction amplitude
		ππ → ππ scattering in S-wave
```
function t(sigma)
     return (g^2)/(SigmaRes- sigma - 1im*rho(sigma)*g^2)
end

function tau(s, sigma)
	(1/3)*Rho(s,sigma)*t(sigma)
end

function BKernel(sigma, s, sigmaPrime)
	((2*s)/(Lambda12(s, sigma)*Lambda12(s, sigmaPrime)))*log((2*s*sigmaPrime +
     2im*s*ϵ - (s + mPi^2 - sigma)*(s + sigmaPrime - mPi^2) -
     Lambda12(s, sigma)*Lambda12(s, sigmaPrime))/(2*s*sigmaPrime +
     2im*s*ϵ - (s + mPi^2 - sigma)*(s + sigmaPrime - mPi^2) +
     Lambda12(s, sigma)*Lambda12(s, sigmaPrime)))
 end
