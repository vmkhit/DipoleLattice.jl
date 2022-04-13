```
  epston(epsilon::Number)
```

```@docs
  Material
```

```@docs
  Interface
```

```@docs
  Layer
```

```@docs
  Stack
```

```@docs
  revert_stack(S::Stack)
```

```@docs
  intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})
```

```@docs
  tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
```

```@docs
  RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
```
```@raw html
<img src="../../Example Pictures/exampleTMM.png" width="60%"/>
```
```@docs
  ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)
```
