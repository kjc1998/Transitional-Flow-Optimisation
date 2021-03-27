


# E_range = 1:0.05:5



function E_Range(E_range)
    blank = zeros(9)
    for i in 1:length(E_range)
        u0 = [0,((E_range[i]-1)^(0.5))/2,((E_range[i]-1)^(0.5))/2,((E_range[i]-1)^(0.5))/2,((E_range[i]-1)^(0.5))/2,0,0,0,0]
        blank = hcat(blank,u0)
    end
    blank = blank[:,2:end]
    return blank
end