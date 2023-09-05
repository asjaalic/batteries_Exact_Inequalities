using DataFrames
using CSV

data = DataFrame(CSV.File("results.csv"))

prices = data[:, :A]
soc = data[:, :SOC]
charge = data[:, :CHARGE]
discharge = data[:, :DISCHARGE]

NSteps = length(prices)

mat = []

for t=1:NSteps

    if charge[t]!=0
        price = prices[t]
        s = soc[t]
        ch = charge[t]
        dis = discharge[t]

        push!(mat, (price,s,ch,dis))

    elseif discharge[t]!=0
        price = prices[t]
        s = soc[t]
        ch = charge[t]
        dis = discharge[t]

        push!(mat, (price,s,ch,dis))

    end

end

new = length(mat)

mat[1][2]

final_matrix = []

for i=1:new-2

    if mat[i+2][4] != 0
        price = mat[i+2][1]+mat[i+1][1]-mat[i][1]
        energy = mat[i+2][4]+mat[i+1][4] + mat[i][3]
        push!(final_matrix, price, energy)
    elseif mat[i+2][3] != 0

    else

    end

end