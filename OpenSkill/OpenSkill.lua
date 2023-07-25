--[[

Tutorial & Documentation: vaschex.github.io/openskill.lua

Version of this module: 1.2.0
Based on philihp/openskill.js commit 291

Ported to Lua & improved by Vaschex

Ported in turn to pure Lua by Totebo 

]]

local module = {}
module.Settings = {
	DefaultModel = "PlackettLuce"
}

--[[
The MIT License (MIT)

Copyright (c) 2020 Philihp Busby <philihp@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
]]

---------------------------------------------------------

local constants = require("OpenSkill.Constants")
local statistics = require("OpenSkill.Statistics")
local util = require("OpenSkill.Util")

local models = {
	["PlackettLuce"] = require("OpenSkill.Models.PlackettLuce"),
	["ThurstoneMosteller"] = require("OpenSkill.Models.ThurstoneMosteller"),
}

function module.Ordinal(rating, options)
	options = options or {}
	return rating.mu - constants.z(options) * rating.sigma
end

function module.Rating(mu, sigma, options)
	options = options or {}
	local copy = {} --copy because i dont want to change the options table
	for k, v in next, options do
		copy[k] = v
	end
	copy.mu = mu
	return {
		mu = mu or constants.mu(options),
		sigma = sigma or constants.sigma(copy)
	}
end

function module.Rate(teams, options)
	options = options or {}
	local model = models[options.model or module.Settings.DefaultModel]
	if options.rank == nil and options.score == nil then
		local result = model(teams, options)
		for i, v in ipairs(result) do
			for i2, v2 in ipairs(v) do
				teams[i][i2].mu = v2[1]
				teams[i][i2].sigma = v2[2]
			end
		end
		return result
	end
	local rank = options.rank
	if not rank then
		rank = {}
		for i, v in ipairs(options.score) do
			rank[i] = -v
		end
	end
	local orderedTeams, tenet = util.unwind(teams, rank)
	local newOptions = {}
	for i, v in next, options do
		newOptions[i] = v
	end
	table.sort(rank)
	newOptions.rank = rank
	local newRatings = model(orderedTeams, newOptions)
	local reorderedTeams = util.unwind(newRatings, tenet)
	for i, v in ipairs(reorderedTeams) do
		for i2, v2 in ipairs(v) do
			teams[i][i2].mu = v2[1]
			teams[i][i2].sigma = v2[2]
		end
	end
	return reorderedTeams
end

function module.WinProbability(teams, options)
	options = options or {}
	local teamRatings = util.teamRating(teams, options)
	local betaSq = constants.betaSq(options)
	local n = #teams
	local denom = (n * (n - 1)) / 2
	local result = {}
	for i, a in ipairs(teamRatings) do
		local prob = 0
		for q, b in ipairs(teamRatings) do
			if i ~= q then
				prob = prob + statistics.phiMajor(
					(a[1] - b[1]) / math.sqrt(n * betaSq + a[2] ^ 2 + b[2] ^ 2))
			end
		end
		table.insert(result, prob / denom)
	end
	return result
end

function module.DrawProbability(teams, options)
	options = options or {}
	local teamRatings = util.teamRating(teams, options)
	local betaSq, beta = constants.betaSq(options), constants.beta(options)
	local n = #teams
	if n == 0 then return
	elseif n == 1 then return 1 end
	
	local div
	if n > 2 then
		div = 1
	else
		div = 2
	end
			
	local denom = (n * (n - 1)) / div
	local playerCount = 0
	for _, v in ipairs(teams) do
		playerCount = playerCount + #v
	end
	local drawMargin = math.sqrt(playerCount) * beta * statistics.phiMajorInverse((1 + 1 / n) / 2)	
	local result = 0
	for i, a in ipairs(teamRatings) do
		for q, b in ipairs(teamRatings) do
			if i ~= q then
				local sigmaBar = math.sqrt(n * betaSq + a[2] ^ 2 + b[2] ^ 2)
				result = result + statistics.phiMajor((drawMargin - a[1] + b[1]) / sigmaBar) -
						statistics.phiMajor((a[1] - b[1] - drawMargin) / sigmaBar)
			end
		end
	end
	return math.abs(result) / denom
end

return module
