using JLD2, Statistics
using LaTeXStrings, Measures, CairoMakie


function cm_to_pt(x)
    x/2.54*72
end

function add_identity_shade(ax; color=:lightgray)
    limits!(ax, ax.finallimits[]) # Do not change the limits anymore
    xlim_ = [ax.xaxis.attributes.limits[]...]
    realxlim = 2 .^ xlim_
    ylim_ = [ax.yaxis.attributes.limits[]...]
    # Find a line of the form y = ax starting from left bottom corner
    y0 = ylim_[1]; x0 = realxlim[1]; b = y0/x0
    bandplot = band!(ax, xlim_, b*realxlim, [y0, y0], color=color)
    translate!(bandplot, 0, 0, -100)
end

function plot_iterations_makie!(f, row, tbl; 
    ylabel="iterations", 
    x = log2.(Ts), 
    active = 1:size(tbl,3), 
    yscale = log10,
    legend = true, 
    showtitles = true,
    showxaxes = true,
    linewidth = 1,
    axiswidth = 0.5linewidth,
    markersize = 6,
    colormap =  Reverse(:viridis)
    )

    syms = [:circle :rect :star5 :diamond :hexagon :utriangle :dtriangle]
    axs = []; lns = []; scs = [];

    for (i, k) = enumerate(active)
        if i==1
            ylabel_ = ylabel
        else 
            ylabel_ = "" 
        end
        if showtitles
            title_ = latexstring("\$\\mathrm{$(fwdAlias[k])}\$")
        else
            title_ = ""
        end
        ax = Axis(f[row,i], 
                  title = title_,
                  xlabel = L"\log_2 T", 
                  ylabel = ylabel_,
                  yscale = yscale,
                  xticks = x, 
                  spinewidth = axiswidth,
                  xtickwidth = axiswidth,
                  ytickwidth = axiswidth,
                  xgridwidth = axiswidth,
                  ygridwidth = axiswidth)
        if i != 1
            hideydecorations!(ax; grid=false)
        end
        if showxaxes == false
            hidexdecorations!(ax; grid=false)
        end
        lns = []; scs = []
        nrows = size(tbl,1)
        for i=1:nrows
            l = Makie.lines!(ax, x, tbl[i,:,k], linewidth=linewidth, color=i, colormap = colormap, colorrange = (1,nrows))
            s = Makie.scatter!(ax, x, tbl[i,:,k], marker=syms[i], markersize=markersize, color=i, colormap = colormap, colorrange = (1,nrows))
            push!(lns, l); push!(scs, s)
        end
        push!(axs, ax)
    end
    linkyaxes!(axs...)
    linkxaxes!(axs...)
    if legend != false
        if legend == true
    Legend(f[1:size(f.layout)[1],size(f.layout)[2]+1], 
      [[l,s] for (l,s) in zip(lns, scs)], 
      map(x -> "$(x-1)", collect(Ns)),
      "N:",
      framevisible = false,
      padding = (2,2,2,2),
      margin = (0,0,0,0),
      patchlabelgap = 2,
      rowgap = -5,
      colgap = 2,
      titlegap = 0,
      groupgap = 0,
    )
        else
            axislegend(axs[end], 
            [[l,s] for (l,s) in zip(lns, scs)], 
            map(x -> "$(x-1)", collect(Ns)),
            "N:",
            framevisible = false,
            nbanks = 2,
            padding = (0,0,0,0),
            margin = (0,0,0,0),
            patchlabelgap = 2,
            rowgap = -7,
            colgap = -2,
            titlegap = -7,
            groupgap = 0
          )
          axs[end].xgridvisible=false
          axs[end].ygridvisible=false
          axs[end].topspinevisible = false
          axs[end].leftspinevisible = false
          axs[end].rightspinevisible = false
          axs[end].bottomspinevisible = false
    end
end
end


function plot_configurations!(f, T, k, a, b; 
    ylabel="iter.", 
    title = "",
    x = log2.(Ts), 
    yscale = log10,
    offset = 0,
    legend = false,
    hideyaxis = false, 
    linear_helper = false,
    linewidth = 1,
    axiswidth = 0.5linewidth,
    markersize = 6,
    colormap = Reverse(:viridis)
    )

    syms = [:circle :rect :star5 :diamond :hexagon :utriangle :dtriangle]
    axs = []
    lns = []; scs = [];
    for i = 1:length(a)
        for j = 1:length(b)
        if j==1
            ylabel_ = latexstring("$ylabel, \$a=$(a[i])\$")
        else 
            ylabel_ = "" 
        end
        if i==1
            if title == ""
                title_ = "b=$(b[j])"
            else
                title_ = latexstring("$title, \$b=$(b[j])\$")
            end
        else
            title_ = ""
        end
        if i==length(a)
            xlabel_ = L"\log_2 T"
        else
            xlabel_ = ""
        end
        ax = Axis(f[i,offset+j], 
                  title = title_,
                  xlabel = xlabel_, 
                  ylabel = ylabel_,
                  yscale = yscale,
                  spinewidth = axiswidth,
                  xtickwidth = axiswidth,
                  ytickwidth = axiswidth,
                  xgridwidth = axiswidth,
                  ygridwidth = axiswidth
                  )
        if j != 1 || hideyaxis
            hideydecorations!(ax; grid=false)
        end
        if i != length(a)
            hidexdecorations!(ax; grid=false)
        end
        lns = []; scs = [];
        tbl = T[i,j]
        nrows = size(tbl,1)
        for m=1:nrows
            l = Makie.lines!(ax, x, tbl[m,:,k], linewidth=linewidth, color=m, colormap = colormap, colorrange = (1,nrows))
            s = Makie.scatter!(ax, x, tbl[m,:,k], marker=syms[m], markersize=markersize, color=m, colormap = colormap, colorrange = (1,nrows))
            push!(lns, l); push!(scs, s)
        end
        push!(axs, ax)
        end
    end
    linkyaxes!(axs...)
    linkxaxes!(axs...)
    if linear_helper
        add_identity_shade.(axs)
    end
    if legend != false
    Legend(f[length(a)+1, 1:6], 
      [[l,s] for (l,s) in zip(lns, scs)], 
      map(x -> "$(x-1)", collect(Ns)),
      "N: ",
      framevisible = false,
      orientation = :horizontal,
      padding = (0,0,0,0),
      margin = (0,0,-5,0),
      patchlabelgap = 2,
      rowgap = -7,
      colgap = 10,
      titlegap = 10,
      groupgap = 0,
      titleposition = :left
    )
    end
end

