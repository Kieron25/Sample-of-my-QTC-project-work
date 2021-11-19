using Yao
using Random
using Plots
using FFTW

Xstr(N::Int) = chain(N, prod([put(N, i=>X) for i=1:N]))
RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.3)) for i=1:N]))
#== or ==#
XstrR(N::Int) = chain(N, repeat(X, 1:N))
Hzz(N::Int) = sum([-1 * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
Uzz(N::Int, Jt::Float64) = time_evolve(Hzz(N), Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true # Can change this to run the alternative XstrR and Rxstr
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> XstrR(N) |> Uzz(N, Jt) |> RXstr(N)
        else
            ψ |> XstrR(N) |> RXstr(N)
        end
        append!(Mz_vec, Yao.expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end
# Number of Qubits and number of cycles is defined here:
q = 8 # N
step = 70 # nsteps

# This bit of code commented out is for plotting a single graph. But I wanted
# to plot multiple graphs seperately so here it is:
#p = 1.2 # Jt - imaginary time evolution?
#title = "DTC plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(p)
#t_vec, Mz_vec = Mz_evolve(q, step, p)
#plt = Plots.plot(t_vec, Mz_vec, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
#Plots.savefig(title*".pdf")

#t_vec = Vector{Float64}()
#Mz_vec = Vector{Float64}()
for i = 1:30
    # 1 after variable names denote they're local variables in the for loop
    title1 = "DTC plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(i/10)
    t_vec1, Mz_vec1 = Mz_evolve(q, step, i/10)
    #append!(t_vec, t_vec1)
    #append!(Mz_vec, Mz_vec1)
    Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
    Plots.title!(title1)
    #Plots.savefig(title1*".png")

    y = fft(Mz_vec1)
    #println(abs.(y), length(y))
    plt = Plots.plot(abs.(y), xlabel = "frequency", xlims = (0, step), ylabel = " \n"*"Fourier components", legend = false)
    #Plots.savefig("Fourier transform for Jt = "*string(i/10)*".png")
    display(plt)
end

# I have no idea what the relevance of this last bit is...
#using YaoExtensions

#c = variational_circuit(6, 6)

#gatecount(c)
#42+72

# Below code is obsolete/ rough code I used to make
# changes to the above code.

#using LsqFit

#ω = [0:10]
#y = fft(Mz_vec)
#println(abs.(y), length(y))
#Plots.plot(abs.(y), xlabel = "frequency", xlims = (0, step), ylabel = " \n"*"Fourier components", legend = false)

#p = [0.9, 1, 0]
#function fitdata(t::Array, p::Array)
#    ft = p[1]*cos(p[2]*t + p[3])
#    return ft
#end

#println(curve_fit(fitdata(t_vec, p), t_vec, Mz_vec, p))
