function KaellenLambda(x,y,z)
	x^2 + y^2 + z^2 - 2*x*y - 2*x*z - 2*y*z
end

function rho(sigma)
	sqrt(Complex(KaellenLambda(sigma,mPi^2,mPi^2)))/(8*pi*sigma)
end

function Lambda12(s, sigma)
	sqrt(Complex(KaellenLambda(s,sigma,mPi^2)))
end

function Rho(s, sigma)
	sqrt(Complex(KaellenLambda(s,sigma,mPi^2)))/(8*pi*s)
end

Kibble(σp,s,σ) = σ*σp*(3*mPi^2+s-σ-σp)-mPi^2*(s-mPi^2)^2
