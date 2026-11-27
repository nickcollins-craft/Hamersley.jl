#=
Tests of boundary conditions: regions and component selections, value prescriptions,
the boundary terms that natural and Robin conditions become, the conflict rule, and
schedules of phases.
=#


"""
    with_conditions(base, conditions; kinematics = base.kinematics, strict = false)

Return `base` with its boundary conditions replaced by `conditions`.
"""
function with_conditions(base, conditions; kinematics=base.kinematics, strict=false)
    problem = Problem(;
        fields=base.fields,
        equations=base.equations,
        boundary_conditions=conditions,
        kinematics=kinematics,
        strict_boundaries=strict,
    )
    return problem
end


"""
    with_schedule(base, schedule)

Return `base` with its boundary conditions replaced by `schedule`.
"""
function with_schedule(base, schedule)
    problem = Problem(;
        fields=base.fields,
        equations=base.equations,
        schedule=schedule,
        kinematics=base.kinematics,
    )
    return problem
end


@testset "Boundary conditions" begin
    elasticity = small_strain_elasticity()
    thermomechanics = finite_strain_thermomechanics()

    @testset "Regions and component selections" begin
        @test Hamersley.as_region(:top) == FacetSet(:top)
        @test Hamersley.region_label(EdgeSet(:top, :front)) ==
            "edge between top and front"
        @test Hamersley.component_selection(nothing) isa Hamersley.AllComponents
        @test Hamersley.component_selection([3, 1, 3]) == Hamersley.SomeComponents((1, 3))

        # Selections overlap without the dimension being known.
        some = Hamersley.SomeComponents
        except = Hamersley.AllComponentsExcept
        @test Hamersley.overlaps(some((1, 2)), some((2,)))
        @test !Hamersley.overlaps(some((1,)), some((2,)))
        @test !Hamersley.overlaps(some((3,)), except((3,)))
        @test Hamersley.overlaps(Hamersley.AllComponents(), some((2,)))

        # Removing components is how inherited conditions are restricted.
        @test Hamersley.remove_components(Hamersley.AllComponents(), some((3,))) ==
            except((3,))
        @test Hamersley.remove_components(except((3,)), some((1,))) == except((1, 3))
        @test Hamersley.remove_components(some((1, 2)), some((1, 2))) === nothing
        @test Hamersley.remove_components(Hamersley.AllComponents(), except((2,))) ==
            some((2,))
    end

    @testset "Value prescriptions" begin
        @test Hamersley.prescription(1.0) == Constant(1.0)
        @test Hamersley.prescription((X, t) -> t) isa Prescribed

        # A table interpolates linearly and holds its end values.
        ramp = Tabulated((0.0, 2.0), (0.0, 10.0))
        @test Hamersley.prescribed_value(ramp, nothing, 0.5) ≈ 2.5
        @test Hamersley.prescribed_value(ramp, nothing, -1.0) == 0.0
        @test Hamersley.prescribed_value(ramp, nothing, 3.0) == 10.0
        vector_ramp = Tabulated((0.0, 1.0), ((0.0, 0.0), (2.0, -4.0)))
        @test all(Hamersley.prescribed_value(vector_ramp, nothing, 0.5) .≈ (1.0, -2.0))
        @test Hamersley.breakpoints(ramp) == (0.0, 2.0)

        # Functions declare their kinks, and may see the current position.
        kinked = Prescribed((X, t) -> abs(t - 1); breakpoints=(1.0,))
        @test Hamersley.breakpoints(kinked) == (1.0,)
        hydrostatic = Prescribed((X, x, t) -> 10.0 - x[2]; current_position=true)
        @test Hamersley.prescribed_value(hydrostatic, (0.0, 0.0), 0.0; x=(0.0, 1.0)) == 9.0

        # Malformed tables and wrongly shaped values are reported.
        @test !isempty(Hamersley.table_issues(Tabulated((1.0, 0.0), (0.0, 1.0))))
        @test !isempty(Hamersley.table_issues(Tabulated((0.0, 1.0), (0.0,))))
        @test isempty(Hamersley.prescription_shape_issues(
            Constant((1.0, 2.0)),
            Hamersley.VectorValue(2),
        ))
        @test !isempty(Hamersley.prescription_shape_issues(
            Constant(1.0),
            Hamersley.VectorValue(nothing),
        ))
    end

    @testset "Checks on single conditions" begin
        # A traction needs a vector, and a scalar field has no components.
        scalar_traction =
            NaturalBC(:u; on=:right, quantity=ReferenceTraction(), value=1.0)
        @test has_issue(verification(with_conditions(elasticity, (scalar_traction,))),
            "needs a tuple of numbers")
        scalar_components = EssentialBC(:T; on=:left, value=0.0, components=1)
        @test has_issue(
            verification(with_conditions(thermomechanics, (scalar_components,))),
            "selects components of a scalar",
        )

        # A selected component may be given as a number, and must have the right length.
        @test passed(verification(with_conditions(
            elasticity,
            (EssentialBC(:u; on=:left, value=0.0, components=2),),
        )))
        @test has_issue(
            verification(with_conditions(
                elasticity,
                (EssentialBC(:u; on=:left, value=(0.0, 0.0), components=2),),
            )),
            "needs a tuple of 1 numbers",
        )

        # Area loads cannot act on points, and only follower loads may see the position.
        point_traction = NaturalBC(
            :u;
            on=VertexSet(:corner),
            quantity=ReferenceTraction(),
            value=(1.0, 0.0),
        )
        @test has_issue(verification(with_conditions(elasticity, (point_traction,))),
            "concentrated forces")
        dead_position = NaturalBC(
            :u;
            on=:right,
            quantity=ReferenceTraction(),
            value=Prescribed((X, x, t) -> (x[1], 0.0); current_position=true),
        )
        @test has_issue(verification(with_conditions(elasticity, (dead_position,))),
            "only a follower")

        # A pinned point of a field with no traces is discrete, so it is not reported.
        u = Field(:u, VectorField(), Lagrange, 2)
        p = Field(:p, ScalarField(), Lagrange, 1; represents=Pressure())
        saddle = Problem(;
            fields=(u, p),
            equations=(
                Equation(
                    :incompressibility,
                    :p,
                    Term(Variation(Value(:p)), Divergence(:u), ReferenceVolume()),
                ),
            ),
            boundary_conditions=(EssentialBC(:p; on=VertexSet(:corner), value=0.0),),
            kinematics=SmallStrain(),
        )
        @test !has_issue(verification(saddle), "does not appear in the weak form")
    end

    @testset "Natural conditions as boundary terms" begin
        # Every natural condition becomes a term of its field's equation.
        @test length(Hamersley.problem_terms(elasticity)) == 2
        @test length(Hamersley.weak_form_terms(elasticity)) == 1

        # A follower pressure makes finite strains nonlinear and couples to the motion.
        pressure = NaturalBC(:u; on=:right, quantity=SpatialPressure(), value=1.0)
        u = Field(:u, VectorField(), Lagrange, 1)
        linear_piola = Constitutive(
            PlaceholderLaw(
                :linear_piola;
                inputs=(F=DeformationGradientTensor(),),
                outputs=(P=FirstPiolaStress(),),
                affine=true,
            ),
            (F=DeformationGradient(:u),),
        )
        momentum = Equation(
            :momentum,
            :u,
            Term(
                Variation(DeformationGradient(:u)),
                output(linear_piola, :P),
                ReferenceVolume(),
            ),
        )
        finite = Problem(;
            fields=(u,),
            equations=(momentum,),
            boundary_conditions=(pressure,),
            kinematics=FiniteStrain(motion=:u),
        )
        @test passed(verification(finite))
        @test !Hamersley.is_affine(finite)
        dead = Problem(;
            fields=(u,),
            equations=(momentum,),
            boundary_conditions=(
                NaturalBC(:u; on=:right, quantity=ReferencePressure(), value=1.0),
            ),
            kinematics=FiniteStrain(motion=:u),
        )
        @test Hamersley.is_affine(dead)

        # A pressure that depends on the current position depends on the motion itself.
        hydrostatic = NaturalBC(
            :u;
            on=:right,
            quantity=SpatialPressure(),
            value=Prescribed((X, x, t) -> 10.0 - x[2]; current_position=true),
        )
        term = Hamersley.condition_term(finite, hydrostatic)
        @test Hamersley.vertex_fields(Hamersley.term_flux(term), finite.kinematics) == [:u]
        @test Hamersley.vertex_fields(Hamersley.term_flux(term), SmallStrain()) == Symbol[]
    end

    @testset "Robin conditions and boundary terms" begin
        convection = RobinBC(:T; on=:outer, coefficient=10.0, reference=293.15)
        robin = with_conditions(thermomechanics, (convection,))
        @test passed(verification(robin))
        statuses = Hamersley.boundary_status(robin)
        @test any(
            status -> status.field === :T && status.status isa Hamersley.TermLoadedTrace,
            statuses,
        )

        # The Robin flux is affine in T, and couples the heat equation to T.
        term = Hamersley.condition_term(robin, convection)
        @test Hamersley.is_affine(Hamersley.term_flux(term), robin.kinematics)
        @test Hamersley.term_dependencies(term, robin.kinematics) == [:T]

        # It may share a face with a natural condition, but not with an essential one.
        heat_flux = NaturalBC(:T; on=:outer, quantity=ReferenceNormalFlux(), value=5.0)
        shared = with_conditions(thermomechanics, (convection, heat_flux))
        @test passed(verification(shared))
        fixed = EssentialBC(:T; on=:outer, value=300.0)
        @test has_issue(
            verification(with_conditions(thermomechanics, (convection, fixed))),
            "claimed by both a Robin condition and an essential condition",
        )

        # A boundary term written in an equation also loads its trace.
        T = Field(:T, ScalarField(), Lagrange, 1; represents=Temperature())
        fourier = Constitutive(
            PlaceholderLaw(
                :fourier;
                inputs=(g=MaterialGradient(ScalarField()),),
                outputs=(q=MaterialFlux(),),
                affine=true,
            ),
            (g=ReferenceGradient(:T),),
        )
        radiation = Constitutive(
            PlaceholderLaw(
                :radiation;
                inputs=(T=Temperature(),),
                outputs=(q=ReferenceNormalFlux(),),
            ),
            (T=Value(:T),),
        )
        heat = Problem(;
            fields=(T,),
            equations=(
                Equation(
                    :heat,
                    :T,
                    Term(
                        Variation(ReferenceGradient(:T)),
                        output(fourier, :q),
                        ReferenceVolume(),
                    ),
                    Term(
                        Variation(Value(:T)),
                        output(radiation, :q),
                        ReferenceBoundary();
                        on=:outer,
                    ),
                ),
            ),
            boundary_conditions=(EssentialBC(:T; on=:inner, value=300.0),),
            kinematics=SmallStrain(),
            strict_boundaries=true,
        )
        @test passed(verification(heat))
        @test occursin("outer: T ValueTrace() — LoadedByTerm", explain(heat))
    end

    @testset "Conflicts" begin
        clamp = EssentialBC(:u; on=:right, value=(0.0, 0.0))
        traction = NaturalBC(:u; on=:right, quantity=ReferenceTraction(), value=(1.0, 0.0))
        @test has_issue(verification(with_conditions(elasticity, (clamp, traction))),
            "claimed by both an essential condition and a natural condition")
        @test has_issue(verification(with_conditions(elasticity, (traction, traction))),
            "claimed by both a natural condition and a natural condition")

        # Different components of one face do not conflict.
        roller = EssentialBC(:u; on=:right, value=0.0, components=1)
        shear = NaturalBC(
            :u;
            on=:right,
            quantity=ReferenceTraction(),
            value=1.0,
            components=2,
        )
        mixed = with_conditions(elasticity, (roller, shear))
        @test passed(verification(mixed))
        @test occursin("right: u ValueTrace() [components 1] — Prescribed", explain(mixed))

        # Two pins of the same point conflict too.
        pin = EssentialBC(:u; on=VertexSet(:corner), value=(0.0, 0.0))
        @test has_issue(verification(with_conditions(elasticity, (pin, pin))),
            "claimed by both an essential condition and an essential condition")
    end

    @testset "Schedules" begin
        # Traction on the right, then the second component held where it stopped.
        clamp = EssentialBC(:u; on=:left, value=(0.0, 0.0))
        traction = NaturalBC(
            :u;
            on=:right,
            quantity=ReferenceTraction(),
            value=Tabulated((0.0, 1.0), ((0.0, 0.0), (1.0, -2.0))),
        )
        hold = EssentialBC(:u; on=:right, value=-0.01, components=2)
        schedule = Schedule(
            Phase(:loading; conditions=(clamp, traction), until=AtTime(1.0)),
            Phase(:holding; conditions=(hold,), until=AtTime(2.0)),
            Phase(:unloading; conditions=(Release(:u; on=:right),)),
        )
        problem = with_schedule(elasticity, schedule)
        @test passed(verification(problem))

        # The traction keeps every component but the one the later phase takes over.
        holding = Hamersley.phase_conditions(problem, :holding)
        kept = only(filter(bc -> bc isa NaturalBC, holding))
        @test Hamersley.bc_components(kept) == Hamersley.AllComponentsExcept((2,))
        @test length(Hamersley.phase_conditions(problem, :unloading)) == 1

        # The report shows each phase.
        text = explain(problem)
        @test occursin("phase holding (until AtTime(2.0)):", text)
        @test occursin("right: u ValueTrace() [all components except 2] — Loaded", text)
        @test occursin("right: u ValueTrace() [components 2] — Prescribed", text)

        # Malformed schedules are reported.
        bad = Schedule(
            Phase(:first; conditions=(clamp,), until=AtTime(2.0)),
            Phase(:first; conditions=(Release(:u; on=:top),), until=AtTime(1.0)),
            Phase(:third),
        )
        report = verification(with_schedule(elasticity, bad))
        @test has_issue(report, "phase first is declared more than once")
        @test has_issue(report, "which is not after 2.0")
        @test has_issue(report, "which it does not inherit")
        endless = Schedule(Phase(:first; conditions=(clamp,)), Phase(:second))
        @test has_issue(verification(with_schedule(elasticity, endless)),
            "is never reached")

        # A conflict in one phase names the phase.
        conflicting = Schedule(
            Phase(:first; conditions=(clamp, traction), until=AtTime(1.0)),
            Phase(:second; conditions=(EssentialBC(:u; on=:left, value=0.0, components=1),
                NaturalBC(:u; on=:left, quantity=ReferenceTraction(), value=1.0,
                    components=1))),
        )
        @test has_issue(verification(with_schedule(elasticity, conflicting)),
            "in phase second, trace ValueTrace() of field u on left")

        # Boundary conditions and a schedule cannot both be given.
        @test_throws ArgumentError Problem(;
            fields=elasticity.fields,
            boundary_conditions=(clamp,),
            schedule=schedule,
            kinematics=SmallStrain(),
        )
        @test_throws ArgumentError Hamersley.phase_conditions(problem, :missing)
    end

    @testset "Linearity per phase" begin
        u = Field(:u, VectorField(), Lagrange, 1)
        linear_piola = Constitutive(
            PlaceholderLaw(
                :linear_piola;
                inputs=(F=DeformationGradientTensor(),),
                outputs=(P=FirstPiolaStress(),),
                affine=true,
            ),
            (F=DeformationGradient(:u),),
        )
        momentum = Equation(
            :momentum,
            :u,
            Term(
                Variation(DeformationGradient(:u)),
                output(linear_piola, :P),
                ReferenceVolume(),
            ),
        )

        # A dead pressure, then a follower pressure on the same face.
        schedule = Schedule(
            Phase(
                :dead;
                conditions=(
                    NaturalBC(:u; on=:right, quantity=ReferencePressure(), value=1.0),
                ),
                until=AtTime(1.0),
            ),
            Phase(
                :follower;
                conditions=(
                    NaturalBC(:u; on=:right, quantity=SpatialPressure(), value=1.0),
                ),
            ),
        )
        problem = Problem(;
            fields=(u,),
            equations=(momentum,),
            schedule=schedule,
            kinematics=FiniteStrain(motion=:u),
        )
        @test passed(verification(problem))
        @test Hamersley.is_affine(problem, :dead)
        @test !Hamersley.is_affine(problem, :follower)
        @test !Hamersley.is_affine(problem)
        @test Hamersley.default_strategy(problem, :dead) isa CachedFactorisationStrategy
        @test Hamersley.default_strategy(problem, :follower) isa NewtonStrategy
    end

    @testset "Regions on grids" begin
        grid = Ferrite.generate_grid(Ferrite.Triangle, (2, 2))
        Ferrite.addvertexset!(grid, "corner", x -> x[1] ≈ -1 && x[2] ≈ -1)
        pinned = with_conditions(
            elasticity,
            (EssentialBC(:u; on=VertexSet(:corner), value=(0.0, 0.0)),),
        )
        @test isempty(Hamersley.grid_issues(grid, pinned))
        missing_vertex = with_conditions(
            elasticity,
            (EssentialBC(:u; on=VertexSet(:nowhere), value=(0.0, 0.0)),),
        )
        @test any(contains("no vertex set named nowhere"),
            Hamersley.grid_issues(grid, missing_vertex))

        # Edges are a 3D region, given as the meeting of two facet sets.
        edge = EssentialBC(:u; on=EdgeSet(:top, :left), value=(0.0, 0.0))
        @test any(contains("needs a 3D grid"),
            Hamersley.grid_issues(grid, with_conditions(elasticity, (edge,))))
        hexahedra = Ferrite.generate_grid(Ferrite.Hexahedron, (1, 1, 1))
        edge_problem = with_conditions(elasticity, (edge,))
        @test isempty(Hamersley.grid_issues(hexahedra, edge_problem))
        named_edge = EssentialBC(:u; on=EdgeSet(:crease), value=(0.0, 0.0))
        @test any(contains("no named edge sets"),
            Hamersley.grid_issues(hexahedra, with_conditions(elasticity, (named_edge,))))
    end
end
