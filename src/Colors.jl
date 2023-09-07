struct ColorsIt
    ls::Vector{Int64}
    max_col::Int64
    rev_dfs::Vector{Int64}
    parents::Vector{Int64}
    has_ci::Bool
    subgraph_ends_rev::Vector{Int64}
    subgraph_ends::Vector{Int64}
    left_siblings::Vector{Int64}
end

function Base.iterate(CI::ColorsIt, col::Vector{Int64}=Int64[])::Union{Nothing,Tuple{NTuple{length(CI.ls),Int64},Vector{Int64}}}

    if isempty(col)
        c::Vector{Int64} = my_minimal_coloring(CI.ls)
        # return ntuple(i -> c[i], length(c)), (init_ColorsIt(CI.ls), c)
        return (c...,), c
    end

    # ((rev_dfs, parents, has_ci, subgraph_ends_rev, subgraph_ends, left_siblings, n), col) = state

    # col::Vector{Int64} = copy(state)

    # has_ci::Bool = get_prop( g, :has_ci )

    # # set up the maximum color, number of vertices, and the parents for the vertices
    # n::Int64 = length( ls )
    # parents::Vector{Int64} = get_prop( g, :parents )

    # # set up the reverse dfs order for the vertices, the end of subgtrees originating from each vertex and 
    # # the left siblings of each vertices
    # rev_dfs::Vector{Int64} = get_prop( g, :reverse_dfs )
    # subgraph_ends::Vector{Int64} = get_prop( g, :ends_of_subgraphs )
    # subgraph_ends_rev::Vector{Int64} = get_prop( g, :ends_of_subgraphs_rev )
    # left_siblings::Vector{Int64} = get_prop( g, :left_siblings )


    #=
        A vertex is colorable if it hasn't reached its maximum possible color.

        The vertex 1 is colorable if either 
        1.) The graph has no central involution and the color is 1 is not the max color; or 
        2.) The graph has central involutiuon and its color is not the second largest color. 
    =#

    # this variable shows how far we need to look in the tree to find the vertex whose 
    # color is modified
    look_end::Int64 = CI.subgraph_ends[1]

    # in colorable we keep the index of the vertex that is to be colored
    colorable::Int64 = 0 #colorable_1 ? 1 : 0

    # we go through the vertices in reverse DFS order
    k::Int64 = 2

    while k <= length(CI.ls)
        v::Int64 = CI.rev_dfs[k]

        # compute end point of the subgraph T_v stemming at vertex v
        v_end::Int64 = CI.subgraph_ends[v]
        v_end_rev::Int64 = CI.subgraph_ends_rev[v]
        # compute the root left sibling of the subgraph T_v
        # (the subgraph stemming from a vertex with the same parent as v, left of v)
        left_s::Int64 = CI.left_siblings[v]

        left_sibling::UnitRange{Int64} = 1:0

        if left_s != -1
            # the left sibling exists
            # we compute the range of indices that correspond to the left sibling
            left_sibling = left_s:(v-1)
        end

        # check if T_v is isomorphic to the left sibling and if it has the same coloring
        # If yes, then the coloring of T_v cannot be increased and hence this part of 
        # the tree can be ignored

        if view(CI.ls, v:v_end) == view(CI.ls, left_sibling) && view(col, v:v_end) == view(col, left_sibling)
            k = findfirst(x -> CI.rev_dfs[x] == left_s, eachindex(CI.ls))
            continue
        end

        # we decide if the color of v can be increased
        if col[v] < CI.max_col - 1
            # we update colorable if necessary 
            if colorable < v
                colorable = v
            end
            #we update look_end 
            look_end = v_end_rev
        elseif col[v] == CI.max_col - 1 && col[CI.parents[v]] != CI.max_col
            if colorable < v
                colorable = v
            end
            look_end = v_end_rev
        end

        if v == look_end && colorable != 0
            # we got the end of a branch that has colorable vertex
            # we need to look no further
            break
        end
        k += 1
    end

    # check if one is colorable 
    colorable_1::Bool = CI.has_ci ? col[1] < CI.max_col - 1 : col[1] < CI.max_col

    # if did not find colorable vertex then return nothing
    if colorable == 0
        if colorable_1
            colorable = 1
        else
            return nothing
        end
    end

    # if colorable == 0 && !colorable_1 
    #     return nothing 
    # elseif colorable == 0 
    #     colorable = 1
    # end

    if colorable != 1 && col[CI.parents[colorable]] == col[colorable] + 1
        # the parent of v already has color col[v] + 1 and so we increase color by two
        col[colorable] += 2
    else
        # else we increase the color by one
        col[colorable] += 1
    end

    # we reset the coloring of each subtree on the right side of v
    # starting from v+1

    j::Int64 = colorable + 1
    while j <= length(CI.ls)
        root_color::Int64 = 0

        # find the end of the subtreee T_j stemming from vertex j
        es::Int64 = CI.subgraph_ends[j]

        # determine the color for the root j of this subtree
        if CI.has_ci && j == 2
            # if has central involution and j is vertex two, then its color
            # is set to the color of vertex 1 plus 1 
            root_color = col[1] + 1
            # else 
            #     # else no constaint on root color
            #     root_color = 0
        end

        # compute the minimal coloring for the subtree T_j
        # with parent color being the color of parent[j]
        col[j:es] = my_minimal_coloring(CI.ls[j:es], parent_color=col[CI.parents[j]], root_color=root_color)

        # the following subtree that needs to be dealt with stems from the end of T_j plus 1 
        j = es + 1
    end

    # return the coloring computed
    return (col...,), col# new_col, ((rev_dfs, parents, has_ci, subgraph_ends_rev, subgraph_ends, left_siblings, n), col)# (g, col)
end

#########AUX FUNCTIONS#######

function init_ColorsIt(ls::Vector{Int64})::Tuple

    n::Int64 = length(ls)

    par::Vector{Int64} = [0 for _ in 1:n]

    foreach(v -> par[v] = findlast(i -> i < v && ls[i] == ls[v] - 1, eachindex(ls)), 2:n)

    stack::Vector{Int64} = [1]
    rev_dfs::Vector{Int64} = []
    v::Int64 = 0
    while length(stack) > 0
        v = pop!(stack)
        push!(rev_dfs, v)
        append!(stack, findall(x -> par[x] == v, 1:n))
    end

    ans = (
        rev_dfs, # reverse_dfs
        par, # parents
        my_has_central_involution(ls), # has_ci
        [end_of_subgraph_rev(ls, rev_dfs, x) for x in 1:n], # ends_of_subgraphs_rev
        [my_end_of_subgraph(ls, x) for x in 1:n], # ends_of_subgraphs
        [my_root_of_left_sibling_subtree(par, x) for x in 1:n] # left_siblings    
    )

    return ans
end

function my_minimal_coloring(ls::Vector{Int64}; parent_color::Int64=0, root_color::Int64=0)::Vector{Int64}#::NTuple{length(ls), Int64}

    # choose the root color

    if root_color == 0
        if parent_color == 0
            root_color = 1
        else
            root_color = parent_color == 1 ? 2 : 1
        end
    end

    if (root_color == 1) == isodd(ls[1])
        col_pair = (2, 1) # = col_even, col_odd
    else
        col_pair = (1, 2)
    end

    ans::Vector{Int64} = [col_pair[(ls[i]%2)+1] for i in eachindex(ls)]
    # for NTuple use ans::NTuple{length(ls),Int64} = col_pair[(ls .% 2) .+ 1] 

    ans[1] = root_color

    return ans
end

function my_root_of_left_sibling_subtree(par::Vector{Int64}, v::Int64)::Int64

    p::Int64 = par[v]
    if p == 0
        return -1
    end
    desc::Vector{Int64} = [x for x in my_children_vertices(par, p) if x < v]

    if isempty(desc)
        return -1
    end

    return maximum(desc)
end

function my_children_vertices(par::Vector{Int64}, v::Int64)::Vector{Int64}

    return findall(i -> par[i] == v, eachindex(par))
end

function my_end_of_subgraph(ls::Vector{Int64}, v::Int64)

    # find the first vertex whose level is not grater then the level of v
    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return findfirst(i -> i == length(ls) || (i >= v && ls[i+1] <= ls[v]), eachindex(ls))
end

function my_end_of_subgraph_rev(ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64)::Int64

    # position of v in r_dfs 
    ps = findfirst(x -> r_dfs[x] == v, 1:length(ls))
    end_ver::Union{Nothing,Int64} = findfirst(k -> ls[k] <= ls[v], view(r_dfs, ps+1:length(r_dfs)))

    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return end_ver === nothing ? r_dfs[end] : r_dfs[ps+end_ver-1]
end

function end_of_subgraph_rev(ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64)::Int64 #to review

    # position of v in r_dfs 
    ps = findfirst(x -> r_dfs[x] == v, 1:length(ls))
    end_ver::Union{Nothing,Int64} = findfirst(k -> ls[k] <= ls[v], r_dfs[ps+1:end])

    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return end_ver === nothing ? r_dfs[end] : r_dfs[ps+end_ver-1]
end

function my_has_central_involution(ls::Vector{Int64})::Bool

    # if the graph is o--o then the answer is yes
    if ls == [1, 2]
        return true
    end

    if iseven(length(ls))
        two = findfirst(i -> i > 2 && ls[i] == 2, eachindex(ls))
        if view(ls, 3:two-1) == 1 .+ view(ls, two:length(ls))
            return true
        end
    end

    return false
end

######Of Iterators#########

function Base.eltype(CI::ColorsIt)
    NTuple{length(CI.ls),Int64}
end
-
### this function counts the isomorphisms of a tree with level sequence ls. Optionally, the tree can be colored with coloration col

function count_iso(ls::Vector{Int64}, col::Tuple{Vararg{Int64}}, m::Tuple{Vararg{Int64}})::Int64

    isempty(m) && return count_iso(ls, col)

    temp_col = Vector{Int64}(undef, length(ls))
    for i in eachindex(temp_col)
        if i in m
            temp_col[i] = typemax(Int64) - i
        else
            temp_col[i] = col[i]
        end
    end

    return count_iso(ls, (temp_col...,))
end

function count_iso(ls::Vector{Int64}, col::Tuple{Vararg{Int64}}=())::Int64

    is_empty::Bool = isempty(col)

    if length(ls) < 3
        if is_empty
            return ls[1] == 1 ? length(ls) : 1 # ls[1] == 1 is equivalent to be a starting-point graph
        else
            return 1
        end
    end

    my_child = findall(i -> i > length(ls) || ls[i] == ls[1] + 1, 1:length(ls)+1)  # find all the children of the root, plus length(ls)+1


    if is_empty
        if ls[1] == 1 # ls[1] == 1 is equivalent to be a starting-point graph
            if iseven(length(ls))  # necessary condition to have the bad involution
                if view(ls, my_child[1]+1:my_child[2]-1) == 1 .+ view(ls, my_child[2]:length(ls))  # check if it has the bad involution
                    a = count_iso(ls[my_child[1]:my_child[2]-1])   # compute the number of colorations of the main subgraph
                    return 2 * (a^2)
                end
            end
        end
    end

    # here we compute the number of colorations of each subgraph
    last_sub::Vector{Int64} = ls[my_child[1]:my_child[2]-1]  # this is the subgraph
    last_col::Tuple{Vararg{Int64}} = col

    if !is_empty
        last_col = col[my_child[1]:my_child[2]-1]
    end

    m::Int64 = 1
    ans::Int64 = 1


    for x in 3:length(my_child)
        if last_sub == view(ls, my_child[x-1]:my_child[x]-1) && (is_empty || last_col == col[my_child[x-1]:my_child[x]-1]) #view(col, my_child[x-1]:my_child[x]-1))
            m += 1
        else
            a = count_iso(last_sub, last_col)   # number of colorations of this subgraph...
            ans *= (a^m) * factorial(m)
            last_sub = ls[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
            # last_col = col[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
            if !is_empty
                last_col = col[my_child[x-1]:my_child[x]-1]
            end
            m = 1
        end
    end

    a = count_iso(last_sub, last_col)  # compute the number of colorations of the last subgraph
    ans *= (a^m) * factorial(m)       # with multiplicity

    return ans
end

#This functions computes the number of colorations. It is not called anywhere in the package, but I decided to keep it because it could be useful in the future

# function n_colorations(ls::Vector{Int64}, n::Int64, magn::Int64)::Int64  # magn is the number of neighbors cones of a cone. It is nc[1]

#     if length(ls) < 3
#         if ls[1] == 1
#             if length(ls) == 1
#                 return n
#             else
#                 return div(n*magn, 2)
#             end
#         else
#             return magn^length(ls)
#         end
#     end

#     my_child = findall(i -> i>length(ls) || ls[i]==ls[1]+ 1, 1:length(ls)+1)  # find all the children of the root, plus length(ls)+1


#     if ls[1] == 1  # ls[1] == 1 is equivalent to be a starting-point graph
#         if iseven(length(ls))  # necessary condition to have the bad involution
#             if view(ls,my_child[1]+1:my_child[2]-1) == 1 .+ view(ls,my_child[2]:length(ls))  # check if it has the bad involution
#                 c = n_colorations(ls[my_child[1]:my_child[2]-1], n, magn)   # compute the number of colorations of the main subgraph
#                 return div(n*(c^2),2*magn)  # return the number of coloration computing only the number of colorations of the main subtree
#             end
#         end
#         ans = n  # in the starting-point graph, the root can assume n values since it has no parent
#     else
#         ans = magn # we are not in the starting-point graph, so the root has a parent to deal with
#     end

#     # here we compute the number of colorations of each subgraph
#     last_sub = ls[my_child[1]:my_child[2]-1]  # this is the subgraph
#     m = 1                                     # this is its multiplicity

#     for x in 3:length(my_child)
#         if last_sub == view(ls, my_child[x-1]:my_child[x]-1)  # we run up to find a different subgraph
#             m += 1
#         else
#             c = n_colorations(last_sub, n, magn)   # number of colorations of this subgraph...
#             ans *= binomial(c+m-1,m)         # ...counted with multiplicity
#             last_sub = ls[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
#             m = 1
#         end
#     end

#     c = n_colorations(last_sub, n, magn)  # compute the number of colorations of the last subgraph
#     ans *= binomial(c+m-1,m)        # with multiplicity

#     return ans
# end
