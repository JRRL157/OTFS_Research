M = 256;            % Modulation order
k = log2(M);       % Bits per symbol
numBits = k*2.5e5; % Total bits to process
sps = 4;           % Samples per symbol (oversampling factor)
filtlen = 10;      % Filter length in symbols
rolloff = 0.25;    % Filter rolloff factor

rng default;                     % Default random number generator
dataIn = randi([0 1],numBits,1); % Generate vector of binary data

constrlen = [5 4];          % Code constraint length
genpoly = [23 35 0; 0 5 13]; % Generator polynomials
tPoly = poly2trellis(constrlen,genpoly);
codeRate = 2/3;

dataEnc = convenc(dataIn,tPoly);

dataSymbolsIn = bit2int(dataEnc,k);
dataMod = qammod(dataSymbolsIn,M);

rrcFilter = rcosdesign(rolloff,filtlen,sps);
txSignal = upfirdn(dataMod,rrcFilter,sps,1);

EbNo = 20;
snr = convertSNR(EbNo,'ebno', samplespersymbol=sps, bitspersymbol=k,CodingRate=codeRate);

rxSignal = awgn(txSignal,snr,'measured');

rxFiltSignal = upfirdn(rxSignal,rrcFilter,1,sps);       % Downsample and filter
rxFiltSignal = rxFiltSignal(filtlen + 1:end - filtlen); % Account for delay

dataSymbOut = qamdemod(rxFiltSignal,M);
codedDataOut = int2bit(dataSymbOut,k);

traceBack = 16;
numCodeWords = floor(length(codedDataOut)*2/3);
dataOut = vitdec(codedDataOut(1:numCodeWords*3/2), tPoly,traceBack,'cont','hard');

decDelay = 2*traceBack;              % Decoder delay, in bits
[numErrors,ber] = biterr(dataIn(1:end - decDelay),dataOut(decDelay + 1:end));       
fprintf('\nThe bit error rate is %5.2e, based on %d errors.\n', ber,numErrors)