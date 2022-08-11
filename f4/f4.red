module f4;

%--------------------------------------------------------------------------------------------------

% Reduce-specific:
sizeofInt32 = 32;

asserted procedure cdrn(x, n);
    if n = 0 then
        x
    else
        cdrn(cdr x, n - 1);    

asserted procedure getvlast(x);
    getv(x, length(x));

%--------------------------------------------------------------------------------------------------

put('f4, 'psopfn, 'f4_groebner);

asserted procedure f4_groebner(u: List): List;
begin scalar x;
    polynomials := reval pop u;
    ring . internalpolys := io_convert_to_internal(polynomials);
    coeffs . exps := internalpolys;
    basis := core_groebner(ring, coeffs, exps);
    return io_convert_to_output(basis)
end;


comment Project structure:

	f4.red is the interface

	gb.red is the implementation
	
	hash.red - hash table with monomials encapsulated
	
	symbolic.red - symbolic prepropcessing step
	
	basis.red - basis and S-pair implementation
	
	modular.red - modular stuff (reduction, reconstructionm correctness checks)
	
	lucky.red - lucky prime numbers manipulations



mirror of `groebner` function, Julia reference:
	
	https://sumiya11.github.io/Groebner.jl/interface/  

Notes on function optional params:
	
	No modular correctness checks
	
	Always interreduce in the end
	
	Ordering is always input
	
	linalg is always exact
		
	certify is always false (for now)
	
	random number generator is always (???)
		
	logging is enabled	
	

Types declarations outline:

	CoeffFF: Integer, CoeffQQ: SQ
	
	ExponentVector: Vector{Integer}
	
	HashKey: Integer, HashValue: Integer
	
	PolyRing: {term order code (lex, deglex, degrevlex): 	INteger, 
					nvars: INteger, 
					explen: Integer}
					
	SPair: {Poly1:Integer, Poly2: Integer, 
				lcm: HashValue, deg: Integer}
				
	Pairset: {Vector{SPair} , load: INteger}
	
	Basis: {ndone:Integer, ntotal: Integer, 
				isredundant: VEctor{Int}, nonredundant: VEctor{Int},
				leaddiv: VEctor{Int}, nlead: Integer };
				
	
		
	
	

% AM f4


% SM f4

function groebner(
            polynomials::Vector{Poly};
            reduced::Bool=true,
            ordering::Symbol=:input,
            certify::Bool=false,
            forsolve::Bool=false,
            linalg::Symbol=:exact,
            rng::Rng=Random.MersenneTwister(42),
            loglevel::Logging.LogLevel=Logging.Warn
            ) where {Poly, Rng<:Random.AbstractRNG}

    #= set the logger =#
    prev_logger = Logging.global_logger(ConsoleLogger(stderr, loglevel))

    #= extract ring information, exponents and coefficients
       from input polynomials =#
    # Copies input, so that polynomials would not be changed itself.
    ring, exps, coeffs = convert_to_internal(polynomials, ordering)

    #= check and set algorithm parameters =#
    metainfo = set_metaparameters(ring, ordering, certify, forsolve, linalg, rng)
    # now ring stores computation ordering
    # metainfo is now a struct to store target ordering

    #= change input ordering if needed =#
    assure_ordering!(ring, exps, coeffs, metainfo)

    #= compute the groebner basis =#
    if ring.ch != 0
        # bexps, bcoeffs = groebner_ff(ring, exps, coeffs, reduced, rng, metainfo)
        # if finite field
        # Always returns UInt coefficients #
        bexps, bcoeffs = groebner_ff(ring, exps, coeffs, reduced, metainfo)
    else
        # if rational coefficients
        # Always returns rational coefficients #
        bexps, bcoeffs = groebner_qq(ring, exps, coeffs, reduced, metainfo)
    end

    # ordering in bexps here matches target ordering in metainfo

    #= revert logger =#
    Logging.global_logger(prev_logger)

    # ring contains ordering of computation, it is the requested ordering
    #= convert result back to representation of input =#
    convert_to_output(ring, polynomials, bexps, bcoeffs, metainfo)
end

endmodule; % end of f4

end;