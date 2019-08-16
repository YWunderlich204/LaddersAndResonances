function balanceplot(z; xy = nothing, cscale::Number = 0.1,
                    title = "")
    zlim = max(abs.(z)...) * cscale
    p = heatmap(z, c=:balance, clim=(-zlim,zlim), title=title)
    # xlab, ylab
    if xy != nothing
        N = size(z)[1]
        Nsteps = 6; step=div(N-1,Nsteps); i0 = div(mod(N-1,Nsteps),2)+1
        xylab = round.(xy, sigdigits=2)[i0:step:N]
        p = heatmap!(z, c=:balance, clim=(-zlim,zlim), title=title,
            xticks=(i0:step:N,xylab), yticks=(i0:step:N,xylab))
    end
    return p
end
