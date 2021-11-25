using Turing

@model function gdemo(x, y)
    s ~ InverseGamma(2, 3)
    m ~ Normal(0, sqrt(s))
    x ~ filldist(Normal(m, sqrt(s)), length(y))
    @show x
    for i in 1:length(y)
        @show x[i]
        y[i] ~ Normal(x[i], sqrt(s))
    end
end

model_gdemo = gdemo([1.0, 0.0], [1.5, 0.0])
c2 = sample(model_gdemo, NUTS(0.65), 100)

result1 = prob"y = [1.5] | chain=c2, model = model_gdemo, x = [1.0]"

@model function simple_model(OBJ)
    a ~ Uniform(0.2,0.8)

    b = a * 10

    OBJ ~ Normal(b, 3.4)
    print(OBJ)
end

init_model = simple_model(5)
c2 = sample(init_model, NUTS(0.65), 100)

result1 = prob"OBJ = 6 | chain=c2, model = init_model, a = 0.4"
