local statistics = {}
local gaussian = require("OpenSkill.Gaussian")

local normal = gaussian.new(0, 1)

function statistics.phiMajor(x)
	return normal:cdf(x)
end

function statistics.phiMajorInverse(x)
	return normal:ppf(x)
end

function statistics.phiMinor(x)
	return normal:pdf(x)
end

function statistics.v(x, t)
	local xt = x - t
	local denom = statistics.phiMajor(xt)
	
	local result
	if denom < 2.2204460492503e-16 then
		result = -xt
	else
		result = statistics.phiMinor(xt) / denom
	end
	
	return result
end

function statistics.w(x, t)
	local xt = x - t
	local denom = statistics.phiMajor(xt)
	if denom < 2.2204460492503e-16 then
		
		local result
		if x < 0 then
			result = 1
		else
			result = 0
		end
		
		return result
	end
	return statistics.v(x, t) * (statistics.v(x, t) + xt)
end

function statistics.vt(x, t)
	local xx = math.abs(x)
	local b = statistics.phiMajor(t - xx) - statistics.phiMajor(-t - xx)
	if b < 1e-5 then
		if x < 0 then return -x - t end
		return -x + t
	end
	local a = statistics.phiMinor(-t - xx) - statistics.phiMinor(t - xx)

	local result
	if x < 0 then
		result = -a
	else
		result = a
	end
	result = result / b
		
	return result
end

function statistics.wt(x, t)
	local xx = math.abs(x)
	local b = statistics.phiMajor(t - xx) - statistics.phiMajor(-t - xx)

	local result
	if b < 2.2204460492503e-16 then
		result = 1
	else
		result = ((t - xx) * statistics.phiMinor(t - xx) +
		(t + xx) * statistics.phiMinor(-t - xx)) / b + statistics.vt(x, t) * statistics.vt(x, t)
	end

	return result
end

return statistics