% Ideally "resamples" x vector from s to u by sinc interpolation
function y = sinc_interp(x,s,u)
    % Interpolates x sampled at "s" uniformly spaced instants
    % Output y is sampled at "u" uniformly spaced instants
    % ("s" for "sampled" and "u" for "upsampled")
    % (consequently, length(x)=length(s))

    % Find the period of the undersampled signal
    T = s(2)-s(1);

    % The entries of this matrix are each u-s permutation.
    % It will be used to generate the sinc transform that will
    % be convolved below with the input signal to do the
    % interpolation.
    %
    % (recall that u(:) will be a column vector regardless
    % of the row-ness of u. So u(:) is a row, and s(:) is a
    % column)
    sincM = repmat( u(:), 1, length(s) ) ...
           - repmat( s(:)', length(u), 1 );

    % * Sinc is the inverse Fourier transform of the boxcar in
    % the frequency domain that was used to filter out the
    % ambiguous copies of the signal generated from sampling.
    % * That sinc, which is now sampled at length(u) instants,
    % is convolved with the input signal becuse the boxcar was
    % multipled with its Fourier transform.
    % So this multiplication (which is a matrix transformation
    % of the input vector x) is an implementation of a
    % convolution.
    % (reshape is used to ensure y has same shape as upsampled u)
    y = reshape( sinc( sincM/T )*x(:) , size(u) );
end