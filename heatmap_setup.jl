using Dates

idx = 1
foldername = "$(string(today()))_run_1"

while ispath("sims/$foldername")
    global foldername = foldername[1:end-length(idx)] * string(idx+1)
    global idx += 1
end
mkpath("sims/$foldername")

open("sims/$foldername/MUpv_phase.csv", "w") do file
    write(file, "Rvv,muvv,mupv,sigmav,C,conv_status\n")
end
open("sims/$foldername/Rvv_phase.csv", "w") do file
    write(file, "Rvv,muvv,mupv,sigmav,C,conv_status\n")
end

println(foldername)