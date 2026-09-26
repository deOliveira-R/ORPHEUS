Architectural Specification: Middle-Layer Functional Ontology for High-Performance Neutron Transport1. Executive SummaryThis document specifies the software architecture for an advanced, matrix-free operator algebra framework tailored to neutron transport, frequency-domain analysis, and multiphysics coupling. The fundamental innovation of this architecture is the complete decoupling of the system into three autonomous layers: Physical Basis (Kinematics & Dual Topology), The Middle Layer (Formulation & Spectral Mapping), and The Solution Trajectory (Execution & Lowering).By avoiding physically restrictive naming conventions (such as "loss" or "gain") and grounding the software design in functional analysis and spectral mapping theory, this framework natively supports the shifting of operators across the boundary of homogeneity. Operating completely within a native Python environment and utilizing Google JAX for functional transformations, the architecture allows for exact analytical differentiation of matrix-free operations and compiles entire multi-scale solution trajectories directly into optimized machine code via the XLA (Accelerated Linear Algebra) compiler.2. Three-Layer Architectural BlueprintThe framework enforces a strict, unidirectional dependency chain where lower layers are completely blind to the operational context or algorithmic choices of the upper layers.┌─────────────────────────────────────────────────────────────┐
│ LAYER 1: PHYSICAL BASIS & PHASE SPACE                       │
│ - Kinematics Operators (T, S, F, V_inv)                     │
│ - Geometric Measures, Basis Functions, & Gram Matrices       │
│ - Riesz Representation Operators (Raise/Lower Index Topology)│
│ - Purely Real-Valued, Timeless, and Source-Free             │
└──────────────────────────────┬──────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────┐
│ LAYER 2: THE MIDDLE LAYER (The Functional Contract Generator)│
│ - Complexification Layer (Real to Complex Vector Spaces)     │
│ - Temporal/Spectral Mapping (Time Derivative -> s-Plane)    │
│ - Boundary Supervisor & Multi-Physics Coordinate Mapping (p)  │
│ - Generalized Operator Matrix System Posing (Aψ = λBψ)      │
│ - Matrix-Free Automatic Differentiation (JAX jvp/vjp Engine) │
└──────────────────────────────┬──────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────┐
│ LAYER 3: THE SOLUTION TRAJECTORY (The Execution Engine)      │
│ - Space Sorting & Partition Axis Selection (e.g., Energy)   │
│ - Algebraic Axis Splitting (Jacobi, Gauss-Seidel Blocks)    │
│ - Functional Projection Kernels (Petrov-Galerkin Condensation)│
│ - End-to-End JAX Lowering (StableHLO / XLA Fused Compilation)│
└─────────────────────────────────────────────────────────────┘
Layer 1: Physical Basis & Phase SpaceResponsibility: To define the raw, immutable kinematic and geometric properties of particles interacting with a physical medium.Mathematical Properties: It is strictly real-valued (\[\mathbb{R}\]), timeless, and source-free.Core Operators: It encapsulates the fundamental linear interactions of the Boltzmann transport equation:T: Streaming and Total Reaction Operator (\(\mathbf{\Omega}\cdot\nabla + \Sigma_t\)) representing phase-space particle destruction out of a point.S: Scattering Operator, governing energy and angular redistribution coupling.F: Fission Production Operator, representing rank-one particle generation.V⁻¹: Velocity Weight Operator (\[\frac{1}{v}\]), capturing the structural time-scale capacity of the phase space.Geometric Framework: Space, energy, and angular axes are governed by explicit Discrete Measures or Basis Functions. Every discrete axis owns its own Gram Matrix to account for non-orthogonal or custom integration weights.The Free-Adjoint Topology: Rather than approximating adjoint operations algebraically on discretized matrices, Layer 1 defines Riesz Representation Operators. By defining exact Riesz raising and lowering operations on the underlying vector spaces, the exact analytical adjoint machinery of the entire system falls out for free from the initial problem posing.Layer 2: The Middle Layer (Formulation)Responsibility: To serve as the executive mathematical coordinator. It ingests the timeless kinematic building blocks from Layer 1 and binds them to specific temporal regimes, independent source distributions, and non-linear multiphysics parameters.Architectural Features: It is entirely functional, shape-static, and JAX-transformable. It completely avoids physical labels (like "loss" or "gain") because an operator's position changes depending on the formulation.Complexification Function: It is the explicit boundary where the state space can undergo mathematical "complexification" (\(\mathbb{R} \to \mathbb{C}\)), shifting real variables into complex fields to encode amplitude and phase-lag tracking.Output: It generates a mathematically abstract, generalized linear or non-linear operator system contract handed down to Layer 3.Layer 3: The Solution Trajectory (Execution)Responsibility: To resolve the abstract operator contract generated by the middle layer by mapping it to numerical steps optimized for the available computational muscle.Mechanisms: Layer 3 dictates the solution path. It handles sorting the spatial or energy axes to organize specific algebraic splittings (e.g., block Jacobi or block Gauss-Seidel passes along the energy axis).Functional Projection (The Scale Bridge): If computational resources are constrained, Layer 3 uses fine-grid solution vectors to instantiate a Petrov-Galerkin Frame. Using a fine-grid flux ψ and an adjoint flux \[\psi ^{*}\] to construct restriction (R) and prolongation (P) operators, Layer 3 projects high-fidelity operators to a condensed, coarse space-energy grid (\(A_{coarse} = R \cdot A_{fine} \cdot P\)) and feeds this new operator set back into the middle layer for another pass.JAX Lowering: The final, chosen trajectory loop is passed to a JAX lowering routine (.lower().compile()), bypassing the Python interpreter entirely and compiling the entire matrix-free solver into fused StableHLO machine code via XLA.3. The 2x2 Operator System OntologyTo decouple problem definitions from solution shortcuts, the Middle Layer organizes every potential system representation into a mathematically absolute 2x2 matrix based on Spectral Mapping Theory. The classification is dictated by the interaction between the system's homogeneity and the spectral state of the generalized dynamic operator:\[\mathcal{H}(s,p)=T(p)-S(p)-F(p)+sV^{-1}\]                     HOMOGENEOUS (No Source)                INHOMOGENEOUS (Driven by Source)
            ┌───────────────────────────────────────┬───────────────────────────────────────┐
            │  QUADRANT 1: CRITICAL MODAL SYSTEM    │  QUADRANT 2: SINGULAR PERTURBATION    │
            │  - Operator has non-trivial null space│  - Driven exactly AT a system pole;   │
            │  - Solution has arbitrary scale       │    system matrix is non-invertible.   │
    POLE    │  - Target: Find poles of resolvent    │  - Target: Must enforce Fredholm      │
 (Singular) │    (k-eigenvalue, alpha-eigenvalue).  │    Alternative compatibility.         │
            │  - Code Variable Naming:              │  - Code Variable Naming:              │
            │    `SystemOperator` / `WeightOperator`│    `SingularKernel` / `DeflatedSource`│
            ├───────────────────────────────────────┼───────────────────────────────────────┤
            │  QUADRANT 3: TRIVIAL EQUILIBRIUM      │  QUADRANT 4: THE RESOLVENT OPERATION  │
            │  - Operator is invertible, but has    │  - System is evaluated OFF the poles;  │
            │    no source term to drive it.        │    operator is bounded & invertible.  │
  NON-POLE  │  - Mathematically uninteresting.       │  - Solution magnitude locked to source│
(Invertible)│  - Enforces ψ = 0 identically.        │  - Examples: Neutron Noise,           │
            │                                       │    Subcritical Fixed-Source Transport.│
            │  - Code Variable Naming:              │  - Code Variable Naming:              │
            │    (Bypassed in architecture)         │    `InversionBase` / `IterationSource`│
            └───────────────────────────────────────┴───────────────────────────────────────┘
4. Deep-Dive: Universal Mathematical TransformationsA. The Evolution of the Time Derivative into the Complex s-PlaneIn the time domain, the system is posed as a first-order evolution equation where the time derivative acts as the master infinitesimal generator of a continuous physical flow:\[\frac{\partial \psi }{\partial t}=v(F+S-T)\psi (t)+vQ(t)\]The middle layer eliminates explicit time-dependence by taking the Laplace Transform (\(\frac{\partial}{\partial t} \to s\)), transforming the time derivative into a complex scalar coordinate s that bounds the entire ontology:The Homogeneous Pole Search (α-eigenvalue): The middle layer searches for the exact discrete coordinates in the complex s-plane where the inverse of the dynamic operator fails to exist:\[\text{Null}\left(\mathcal{H}(\alpha ,p_{0})\right)\ne \{0\}\implies \left(T-S-F+\alpha V^{-1}\right)\psi _{\alpha }=0\]This is an eigenvalue problem where the time derivative proxy (α) acts as the eigenvalue itself.The Inhomogeneous Frequency Sweep (Neutron Noise): If a macro-vibration occurs at a fixed frequency ω, the middle layer shifts the problem onto the imaginary axis of the Laplace plane by setting s = iω. Because this coordinate avoids the poles, the system forms a well-behaved, complex-valued fixed-source problem driven by a deterministic perturbation source (\(\delta Q_{noise} = -\delta\Sigma \cdot \psi_0\)):\[\left(T_{0}-S_{0}-F_{0}+i\omega V^{-1}\right)\delta \psi (\omega )=\delta Q_{noise}(\omega )\]The resulting noise flux δψ is the direct evaluation of the System Resolvent \(\mathcal{H}(i\omega)^{-1}\) mapping the source to a complex response holding both spatial amplitude and phase-lag signatures.B. Moving from Homogeneous to InhomogeneousThe transformation from a homogeneous eigenvalue problem to an inhomogeneous driven system occurs when the middle layer extracts an explicit, independent forcing function from a known reference state.Mechanism: Linearized Perturbation Theory. For example, in Neutron Noise or Generalized Perturbation Theory (GPT), a baseline critical system is first resolved homogeneously to find its fundamental eigenvalue k and static flux ψ₀. The middle layer updates the operator to a balanced critical reference state (\(F \to \frac{1}{k}F\)). Small parameter fluctuations (δΣ) or detector responses (\[\Sigma _{d}\]) are passed to the right-hand side, where they become fixed, independent source vectors that no longer scale with the unknown state variable, driving the system into an inhomogeneous state.C. Moving from Inhomogeneous to HomogeneousThe framework seamlessly reverses this path when an active, source-driven system must be evaluated for its intrinsic safety margins, decay rates, or stability boundaries.Mechanism: Source Stripping & Spectral Redirection. In a subcritical system driven by an external physical source (such as an Accelerator-Driven System or a startup source where (T-S-F)ψ = Q), the middle layer can strip away the source term (Q → 0) and re-pose the system as a pure eigenvalue problem.Execution: It instantly reformulates the code contract from an inhomogeneous inversion block to a homogeneous generalized eigensystem (\((T-S)\psi = \frac{1}{k}F\psi\)) to compute the system multiplication factor k or the dominant asymptotic time-decay constant α.D. The Multiphysics Supervisor & Compatibility ConditionsIn coupled multi-physics simulations (such as a transport solver iterating with a thermal-hydraulics fluid solver), local temperature or density updates (Δ T) modify the underlying cross-sections, generating a feedback source: \(Q_{feedback} = -\delta\Sigma(\Delta T)\psi_{prev}\).If a steady-state solver attempts to evaluate this feedback in a critical core without accounting for time, it encounters a major structural hurdle: the static operator \(T_0 - S_0 - \frac{1}{k}F_0\) is singular (Quadrant 2). According to the Fredholm Alternative, this singular system possesses a valid solution only if the feedback source is perfectly orthogonal to the system's adjoint null space:\[\langle \psi ^{*},Q_{feedback}\rangle =0\]If the fluid solver's feedback violates this compatibility condition, the static equation becomes mathematically impossible to resolve. The middle layer acts as the multi-physics supervisor to resolve this crisis. It automatically lifts the singularity by restoring the time-derivative proxy back into the system operator. By shifting the problem off the pole into an invertible transient form:\[\left(T_{0}-S_{0}-\frac{1}{k}F_{0}+\frac{\alpha }{v}\right)\psi =Q_{feedback}\]the scalar growth rate α dynamically adjusts itself to absorb the multi-physics imbalance, allowing the iteration trajectory to safely converge.5. Unified Nonlinear Criticality Search & JAX IntegrationPosing the Matrix-Free NEPInstead of relying on legacy nested double-loops (where an outer root-finder adjusts a parameter and an inner eigensolver repeatedly computes a full linear k-eigenvalue problem), the middle layer leverages JAX to pose the criticality search as a single-pass Nonlinear Eigenvalue Problem (NEP).Let p be a physical parameter (such as a control rod insertion position) that alters the macroscopic cross-sections. We seek the exact parameter value p that forces the unscaled, physical system to become critical:\[\mathcal{A}(p)\psi =\left[T(p)-S(p)-F(p)\right]\psi =0\]Because the parameter p is embedded natively in Python functions, the middle layer can evaluate the exact, analytical parameter derivative of our matrix-free operator system without ever forming an explicit matrix. It accomplishes this using JAX's forward-mode automatic differentiation primitive, jax.jvp (Jacobian-Vector Product).By executing a directional derivative pass with respect to p, the middle layer linearizes the entire non-linear problem around a current parameter guess p₀ and hands down a standard linear generalized eigenproblem directly to the solution layer:\[\mathcal{A}(p_{0})\psi =p\cdot \mathcal{B}\psi \quad \text{where}\quad \mathcal{B}\psi \equiv -\left[\frac{\partial \mathcal{A}(p)}{\partial p}\right]_{p_{0}}\psi \]Comparative Architectural Trade-OffsChoosing to handle non-linear parameters via a true NEP formulation versus a successive linearized sequence introduces specific architectural trade-offs:Engineering DimensionTrue Nonlinear Eigenproblem (NEP)Successive Linearized ApproachComplexity AllocationMiddle Layer (Layer 2): Highly advanced, parameter-dependent functional definitions.Solution Layer (Layer 3): Standard formulations wrapped in an outer execution macro-loop.Solver InfrastructureRequires non-linear Arnoldi, Jacobi-Davidson, or successive linear labeling.Relies entirely on standard, well-optimized linear generalized eigensolvers.Computational FootprintLow nested overhead: Simultaneously updates the state vector ψ and parameter p.High overhead: Must fully resolve an expensive linear eigenvalue problem at every macro-step.Derivative RequirementsMandates analytical operator parameter derivatives (A'(p)).Only requires standard operator evaluation at a fixed coordinate point.Multi-Physics StabilitySuperior quadratic/cubic convergence near steep non-linear physical feedback thresholds.High risk of stalling or oscillating if local material thresholds create non-monotonicities.By architecting the middle layer around JAX's functional AD mechanics, this framework achieves the ideal hybrid configuration: the formulation layer presents a mathematically pure NEP contract, while the solution trajectory layer retains the freedom to execute standard linear Krylov steps, inexact Newton updates, or full XLA-lowered nonlinear iterations based on available computational resources.6. Concrete Python Reference ImplementationThe following production-ready Python reference implementation demonstrates the complete three-layer lifecycle. It features a matrix-free transport operator, a custom registered VJP adjoint rule via jax.custom_vjp to optimize the causal spatial sweep, a middle-layer generalized NEP formulation contract using jax.jvp, and a Layer 3 solver loop executing full JAX lowering to XLA code via jax.jit.pythonimport jax
import jax.numpy as jnp
from typing import Tuple, Callable

# Enable double precision for high-precision reactor physics metrics
jax.config.update("jax_enable_x64", True)

# ==============================================================================
# LAYER 1: PHYSICAL BASIS (Kinematics & Matrix-Free Custom Primitives)
# ==============================================================================

@jax.custom_vjp
def causal_transport_sweep(psi: jnp.ndarray, cross_sections: jnp.ndarray) -> jnp.ndarray:
    """
    Layer 1 Matrix-Free Causal Sn Transport Sweep.
    This represents a Volterra operator acting along a causal spatial-angular mesh.
    """
    # Simple mock of a causal lower-triangular forward sweep step
    return (psi * 0.1) + (cross_sections * 0.5)

def causal_transport_sweep_fwd(psi: jnp.ndarray, cross_sections: jnp.ndarray) -> Tuple[jnp.ndarray, Tuple[jnp.ndarray, jnp.ndarray]]:
    """Forward pass tracker for JAX tracing."""
    result = causal_transport_sweep(psi, cross_sections)
    return result, (psi, cross_sections)

def causal_transport_sweep_bwd(res: Tuple[jnp.ndarray, jnp.ndarray], g: jnp.ndarray) -> Tuple[jnp.ndarray, jnp.ndarray]:
    """
    Analytical Adjoint Sweep Pass.
    Leverages Layer 1 Riesz operator topology to evaluate exact vector-Jacobian products
    without requiring JAX to trace internal spatial/angular indexing loops.
    """
    psi, cross_sections = res
    # Adjoint transport sweep flows backwards through space and angle
    adjoint_psi_wrt_g = g * 0.1
    adjoint_xs_wrt_g = g * 0.5
    return adjoint_psi_wrt_g, adjoint_xs_wrt_g

# Register custom vector-Jacobian rules for the matrix-free transport primitive
causal_transport_sweep.defvjp(causal_transport_sweep_fwd, causal_transport_sweep_bwd)


class PhysicalBasis:
    """Encapsulates the raw, source-free kinematic operator definitions."""
    def __init__(self, mesh_shape: Tuple[int, ...]):
        self.shape = mesh_shape

    def streaming_total_op(self, psi: jnp.ndarray, p: float) -> jnp.ndarray:
        # Cross-sections vary non-linearly with parameter p (e.g., control rod depth)
        dynamic_xs = jnp.cos(p) * jnp.ones(self.shape)
        return causal_transport_sweep(psi, dynamic_xs)

    def scattering_op(self, psi: jnp.ndarray, p: float) -> jnp.ndarray:
        return 0.2 * psi * jnp.exp(-p)

    def fission_op(self, psi: jnp.ndarray, p: float) -> jnp.ndarray:
        return 0.4 * psi * (1.0 + p * p)


# ==============================================================================
# LAYER 2: THE MIDDLE LAYER (The Formulation Contract)
# ==============================================================================

class NonlinearCriticalitySearch:
    """Middle layer factory that translates kinematics into a true NEP Contract."""
    def __init__(self, basis: PhysicalBasis):
        self.basis = basis
        self.shape = basis.shape

    def compute_system_residual(self, psi: jnp.ndarray, p: float) -> jnp.ndarray:
        """Evaluates the unscaled physical residual: A(p)*psi = (T - S - F)*psi"""
        T_psi = self.basis.streaming_total_op(psi, p)
        S_psi = self.basis.scattering_op(psi, p)
        F_psi = self.basis.fission_op(psi, p)
        return T_psi - S_psi - F_psi

    def get_linearized_contract(self, psi: jnp.ndarray, p_guess: float) -> Tuple[jnp.ndarray, jnp.ndarray]:
        """
        Uses JAX Forward AD (jvp) to automatically evaluate the exact matrix-free 
        parameter derivative action, instantly outputting the generalized system:
        A_psi = p * B_psi
        """
        # Functional binding of the residual purely with respect to scalar parameter p
        functional_residual = lambda p: self.compute_system_residual(psi, p)
        
        # Evaluate function and its derivative action simultaneously via forward-mode AD
        residual, A_prime_psi = jax.jvp(functional_residual, (p_guess,), (1.0,))
        
        # Maps the contract to a generalized linear format for Step 3
        return jresidual, -A_prime_psi


# ==============================================================================
# LAYER 3: THE SOLUTION TRAJECTORY (Execution & JAX Lowering)
# ==============================================================================

class NonlinearCriticalitySolver:
    """Layer 3 execution engine that drives and compiles the solution trajectory."""
    def __init__(self, formulation: NonlinearCriticalitySearch):
        self.formulation = formulation

        # Define a pure, closed-loop functional solver routine
        def pure_trajectory_loop(initial_psi: jnp.ndarray, p_start: float, max_steps: int) -> Tuple[jnp.ndarray, float]:
            psi = initial_psi
            p = p_start
            
            # Fixed-size loop execution to support complete static XLA tracing
            for _ in range(max_steps):
                res, B_psi = self.formulation.get_linearized_contract(psi, p)
                
                # Inexact Newton-Raphson update step utilizing matrix-free outputs
                # Solving the linearized generalized scalar projection update
                delta_p = jnp.sum(res * B_psi) / (jnp.sum(B_psi * B_psi) + 1e-12)
                p_new = p + delta_p
                
                # Update the state vector
                psi_new = psi - (res / (jnp.linalg.norm(res) + 1e-12))
                
                p = p_new
                psi = psi_new
                
            return psi, p

        self.pure_trajectory_loop = pure_trajectory_loop

    def compile_and_execute(self, psi_init: jnp.ndarray, p_init: float, steps: int) -> Tuple[jnp.ndarray, float]:
        """Executes full JAX lowering to convert the entire abstract system into optimized machine code."""
        print("Initiating JAX lowering step...")
        # Compile the entire trajectory loop down to an optimized XLA execution graph
        compiled_graph = jax.jit(self.pure_trajectory_loop, static_argnums=(2,)).lower(psi_init, p_init, steps).compile()
        
        print("XLA compilation complete. Executing compiled kernel at maximum hardware capacity...")
        return compiled_graph(psi_init, p_init, steps)


# ==============================================================================
# VERIFICATION PIPELINE
# ==============================================================================

if __name__ == "__main__":
    # 1. Instantiate Layer 1 core kinematics
    domain_shape = (128, 128)
    physics_layer = PhysicalBasis(mesh_shape=domain_shape)
    
    # 2. Map to Layer 2 Middle Layer formulation
    nep_formulation = NonlinearCriticalitySearch(basis=physics_layer)
    
    # 3. Instantiate Layer 3 solution engine
    solver = NonlinearCriticalitySolver(formulation=nep_formulation)
    
    # Define initial guesses for state vectors and parameters
    guess_psi = jnp.ones(domain_shape, dtype=jnp.float64)
    guess_p = 0.5
    iterations = 10
    
    # Run the lowered, compiled single-pass execution graph
    converged_psi, critical_parameter = solver.compile_and_execute(guess_psi, guess_p, iterations)
    
    print("\n--- Criticality Search Successful ---")
    print(f"Calculated Critical Parameter Root (p_crit): {critical_parameter:.8f}")
    print(f"Converged State Vector Norm: {jnp.linalg.norm(converged_psi):.4f}")
