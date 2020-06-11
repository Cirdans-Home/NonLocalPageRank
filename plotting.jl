using CSV, DataFrames, PyCall,  Seaborn, LaTeXStrings, Pandas

pygui(:tk)

close("all")



figure()
df = CSV.read("./results/results_cv_linkpred_june11.csv")

res = []
for row in eachrow(df)
    global res;
    r1 = (dataset=row[:dataset], num_nodes=row[:num_nodes], num_edges=row[:num_edges], density=row[:density], trial=row[:trial], alpha=row[:alpha], c=row[:c_nonlocal], accuracy=row[:accuracy_nonlocal], method="nonlocal")
    r2 = (dataset=row[:dataset], num_nodes=row[:num_nodes], num_edges=row[:num_edges], density=row[:density], trial=row[:trial], alpha=0, c=row[:c_local], accuracy=row[:accuracy_local], method="local")
    push!(res, r1,r2)
end
newdf = Pandas.DataFrame(res)


Seaborn.boxplot(x="dataset", y="accuracy", hue="method", data=newdf)
