#=
Shared test problems.

Small formulations used across several test files: small-strain linear elasticity, and
finite-strain thermo-mechanics with a spatial (Fourier) heat flux that is pulled back.
=#


"""
    small_strain_elasticity(; linear = nothing)

Small-strain linear elasticity, ``∫ σ : δε \\, \\mathrm{d}v_0``, clamped on `left` and
loaded by a traction on `right`.
"""
function small_strain_elasticity(; linear=nothing)
    u = Field(:u, VectorField(), Lagrange, 1)
    elasticity = Constitutive(
        IsotropicLinearElasticity(first_lame_parameter=1.0, shear_modulus=1.0),
        (ε=InfinitesimalStrain(:u),),
    )
    momentum = Equation(
        :momentum,
        :u,
        Term(Variation(InfinitesimalStrain(:u)), output(elasticity, :σ), ReferenceVolume()),
    )
    problem = Problem(;
        fields=(u,),
        equations=(momentum,),
        boundary_conditions=(
            EssentialBC(:u; on=:left, value=PlaceholderValue(:u_bar)),
            NaturalBC(:u; on=:right, quantity=ReferenceTraction()),
        ),
        kinematics=SmallStrain(),
        linear=linear,
    )
    return problem
end


"""
    finite_strain_thermomechanics(; traction = ReferenceTraction())

Finite-strain thermo-mechanics: ``∫ P(F, T) : δF \\, \\mathrm{d}v_0`` for momentum, and
``∫ (J F^{-1} q) \\cdot ∇_0 δT \\, \\mathrm{d}v_0`` for heat, with ``q`` given by a spatial
Fourier law in ``∇_x T = F^{-T} ∇_0 T``.
"""
function finite_strain_thermomechanics(; traction=ReferenceTraction())
    u = Field(:u, VectorField(), Lagrange, 2)
    T = Field(:T, ScalarField(), Lagrange, 1; represents=Temperature())
    F = DeformationGradient(:u)

    # The mechanical law depends on both F and T, so momentum couples to temperature.
    mechanics = Constitutive(
        PlaceholderLaw(
            :thermo_hyperelastic;
            inputs=(F=DeformationGradientTensor(), T=Temperature()),
            outputs=(P=FirstPiolaStress(),),
        ),
        (F=F, T=Value(:T)),
    )

    # The Fourier law is spatial, so its input is pushed forward and its output pulled back.
    fourier = Constitutive(
        PlaceholderLaw(
            :fourier;
            inputs=(g=SpatialGradient(ScalarField()),),
            outputs=(q=SpatialFlux(),),
        ),
        (g=PushForward(ReferenceGradient(:T)),),
    )

    problem = Problem(;
        fields=(u, T),
        equations=(
            Equation(
                :momentum,
                :u,
                Term(
                    Variation(F),
                    output(mechanics, :P),
                    ReferenceVolume(),
                ),
            ),
            Equation(
                :heat,
                :T,
                Term(
                    Variation(ReferenceGradient(:T)),
                    PullBack(output(fourier, :q)),
                    ReferenceVolume(),
                ),
            ),
        ),
        boundary_conditions=(
            NaturalBC(:u; on=:right, quantity=traction),
            EssentialBC(:T; on=:left, value=PlaceholderValue(:T_bar)),
        ),
        kinematics=FiniteStrain(motion=:u),
    )
    return problem
end
