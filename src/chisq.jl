function chisq_goodness(x::AbstractVector{T}) where {T <: Real}
    # figure out number of bins, which depends on length(x)
    nbins = round(Int, ceil(2 * (length(x)^(2/5))))

    # create a `range` object
    bin_range = range(-3, 3, length = nbins + 1)

    # use StatsBase `fit` function to count number of observations per bin
    bincounts = fit(Histogram, x, bin_range).weights

    # use `pdf` function to find expected number of observations per bin
    bin_centers = bin_range[1:end-1] .+ bin_range.step.hi / 2
    expected_pdf = pdf.(Normal(0, 1), bin_centers)
    bin_expected = floor.(Int, expected_pdf .* sum(bincounts) / sum(expected_pdf))

    # calcuate pearson stats (X2)
    𝜒² = pearson_chi(bincounts, bin_expected)

    # calculate X2/DF and return it
    return 𝜒² / nbins
end

function pearson_chi(
    observed::AbstractVector{T},
    expected::AbstractVector{T}
) where {T <: Integer}
    χ² = 0
    for i = 1:length(observed)
        χ² += (observed[i] - expected[i])^2 / expected[i]
    end
    return χ²
end
