clc; 
clear all; 
close all;

%% Ask user to select an image
 [filename, pathname] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.tif;*.webp', 'Image Files'}, 'Select an Image');
 if isequal(filename,0) || isequal(pathname,0)
     disp('No image selected. Exiting.');
     return;
 else
     selected_file = fullfile(pathname, filename);
     disp(['User selected ', selected_file])
 end

%% Load and process the image
  img_path = selected_file;
  original_image = imread(img_path);

% Convert to grayscale if needed
  if size(original_image, 3) == 3
      original_image = rgb2gray(original_image);
  end
    
% Display original image
  figure('Name', 'Original Image', 'NumberTitle', 'on');
  imshow(original_image);
  title('Original Grayscale Image');

% Get image dimensions
  [rows, cols] = size(original_image);
  fprintf('Image dimensions: %d x %d\n', rows, cols);
  fprintf('Total pixels: %d\n', rows * cols);

%% Load and process the image
  img_vector = original_image(:);    % Stack columns into a single column
   
% Calculate probabilities of each symbol 
  symbols = unique(img_vector); % Remove repeated symbols 
  probabilities = zeros(size(symbols)); % Initialize probability array 
  
  for i = 1:length(symbols) 
      probabilities(i) = sum(img_vector == symbols(i)) / length(img_vector); 
  end 

% Calculate information for each symbol 
  info = -log2(probabilities); 
  info_to_ent = probabilities .* info; 

% Calculate entropy 
  entropy = sum(info_to_ent(isfinite(info_to_ent))); 

% Display information and entropy
  fprintf('%-12s | %-12s | %-16s\n','Gray level(Symbols)','Probability','Information gained');
  fprintf('---------------------------------------------------------------\n');

  for i = 1:length(symbols)
      fprintf('%-19d | %-12.5f | %-16.4f\n', symbols(i), probabilities(i), info(i));
  end

  fprintf('---------------------------------------------------------------\n');
  fprintf('Entropy: %.6f bits/symbol\n', entropy);

%% Huffman Coding 

% As we deal with bits so no need for dummy symbols : k = ceil(q-r/r-1) >>>
% k = ceil(q-2/2-1) = ceil(q-2/1) = q-2 >>>> so there is no fraction

num_of_symbols = length(symbols); 
codeword = cell(1, num_of_symbols); % Empty cell array for Huffman codewords 
x = zeros(num_of_symbols, num_of_symbols); % Matrix for Huffman algorithm

% Temporary probability array to keep original values intact 
temp_p = probabilities; 

% Huffman algorithm 
for a = 1:num_of_symbols-1 
[~, idx] = sort(temp_p); % Sort probabilities 
% Combine the smallest two probabilities 
temp_p(idx(2)) = temp_p(idx(1)) + temp_p(idx(2));  
x(idx(1), a) = 11; % Assign value to the first smallest p 
x(idx(2), a) = 10; % Assign value to the second smallest p 
temp_p(idx(1)) = nan; % Mark the combined node 
end 

% Generate Huffman codewords based on the matrix 
for a = num_of_symbols-1:-1:1 
code_1 = find(x(:, a) == 11); %return index of 11
code_0 = find(x(:, a) == 10); %return index of 10
codeword{code_1} = [codeword{code_0}, '1']; 
codeword{code_0} = [codeword{code_0}, '0']; 
end 

% Display Huffman dictionary with ASCII code 
fprintf('---------------------------------------------------------------\n'); 
disp('Huffman Coding: Symbol | Codeword'); 
for i = 1:length(symbols) 
fprintf('%d | %s\n', symbols(i), codeword{i}); 
end 
% Encode the data using Huffman codes 
encoded_huffman = strings(1, length(img_vector)); 
for i = 1:length(img_vector) 
for j = 1:length(symbols) 
if img_vector(i) == symbols(j) 
encoded_huffman(i) = codeword{j}; 
end 
end 
end 

% Transmitted Binary Code (Encoded Message - Huffman) 
str = strjoin(encoded_huffman); 
str = strrep(str, ' ', ''); 
disp('---------------------------------'); 
disp('Transmitted Binary Code (Encoded Message - Huffman):'); 
disp(str); 

% Decode the message using Huffman codes 

% Build a reverse dictionary (codeword → symbol) using containers.Map
% This allows very fast decoding by directly mapping each codeword to its symbol.
  revDict_Huffman = containers.Map();

  for k = 1:length(symbols)
    revDict_Huffman(codeword{k}) = symbols(k);
  end

  decoded_Huffman = zeros(1, length(img_vector)); % preallocate
  writeIndex = 1;


  for i = 1:length(encoded_huffman)
    cw = encoded_huffman(i);      % one codeword is complete already
    if isKey(revDict_Huffman, cw)
      symbol = revDict_Huffman(cw);       % Direct lookup
    else
      error('Unknown codeword encountered: %s', cw);
    end         
    decoded_Huffman(writeIndex) = symbol;
    writeIndex = writeIndex + 1;
  end


  disp('---------------------------------'); 
  disp('Decoded Data (Decoded Message - Huffman):'); 
  %disp(decoded_Huffman);

% ======== Display Image ========
  huffman_img = reshape(uint8(decoded_Huffman), size(original_image));
  % Grayscale images in MATLAB use pixel values in the form: datatype: uint8 || range: 0 to 255
  % reshape changes the shape of the array without changing the data order.
  figure('Name', 'Huffman Decoded Image', 'NumberTitle', 'on');
  imshow(huffman_img);
  title('Huffman Decoded Grayscale Image');

%% Shannon Coding using alpha method 
% Sort the symbols and probabilities based on probabilities 
  [sorted_probabilities, idx] = sort(probabilities, 'descend'); 
  sorted_symbols = symbols(idx); % Rearranges the symbols to match their sorted probabilities.

% Calculate alpha values 
  alpha = zeros(1, length(sorted_probabilities)); 
  alpha(1) = 0; % First alpha value is 0 
  for j = 2:length(sorted_probabilities) 
    alpha(j) = alpha(j-1) + sorted_probabilities(j-1); 
  end 

% Calculate code lengths 
  codeLengths = ceil(-log2(sorted_probabilities));

% Compute Shannon codes for each symbol 
  shannon_codeword = cell(1, length(sorted_symbols)); 
  for i = 1:length(sorted_symbols) 
    int = alpha(i); 
    code = []; 
    for j = 1:codeLengths(i) 
      frac = int * 2; 
      c = floor(frac); 
      frac = frac - c; 
      code = [code c]; 
      int = frac; 
    end 
  shannon_codeword{i} = num2str(code); % Convert code array to string 
  end 

% Display Shannon dictionary with alpha values 
  disp('----------------------------------------------'); 
  disp('Shannon Coding');
  fprintf('%-8s | %-12s | %-16s\n','Symbols','Alpha','Codeword');
  for i = 1:length(sorted_symbols) 
    fprintf('%-8d | %.10f | %s\n', sorted_symbols(i), alpha(i), shannon_codeword{i}); 
  end 

% Encode the data using Shannon codes 
  encoded_shannon = strings(1, length(img_vector)); 
  for i = 1:length(img_vector) 
    for j = 1:length(sorted_symbols) 
      if img_vector(i) == sorted_symbols(j) 
        encoded_shannon(i) = shannon_codeword{j}; 
      end 
    end 
  end 

% Transmitted Binary Code (Encoded Message - Shannon) 
  str_shannon = strjoin(encoded_shannon); 
  str_shannon = strrep(str_shannon, ' ', ''); 
  disp('---------------------------------'); 
  disp('Transmitted Binary Code (Encoded Message - Shannon):'); 
  disp(str_shannon); 

% Decode the message using Shannon codes

% Build a reverse dictionary (codeword → symbol) using containers.Map
% This allows very fast decoding by directly mapping each codeword to its symbol.
  revDict = containers.Map();

  for k = 1:length(sorted_symbols)
    revDict(shannon_codeword{k}) = sorted_symbols(k);
  end

  decoded_shannon = zeros(1, length(img_vector)); % preallocate
  writeIndex = 1;


  for i = 1:length(encoded_shannon)
    cw = encoded_shannon(i);      % one codeword is complete already
    if isKey(revDict, cw)
      symbol = revDict(cw);       % Direct lookup
    else
      error('Unknown codeword encountered: %s', cw);
    end         
    decoded_shannon(writeIndex) = symbol;
    writeIndex = writeIndex + 1;
  end

  disp('---------------------------------'); 
  disp('Decoded Data (Decoded Message - Shannon):'); 
  %disp(decoded_shannon);

% ======== Display Image ========
  shannon_img = reshape(uint8(decoded_shannon), size(original_image));
  figure('Name', 'shannon Decoded Image', 'NumberTitle', 'on');
  imshow(shannon_img);
  title('Shannon Decoded Grayscale Image');

%% Display Shannon's code average length and efficiency

% Calculate the average code length for Shannon 
  average_shannon_length = sum(sorted_probabilities .* codeLengths);
  disp('---------------------------------'); 
  fprintf('Code average length (Shannon): %.4f bits/symbol\n', average_shannon_length); 
% Calculate Shannon efficiency (should not exceed 100%) 
  efficiency_shannon = (entropy / average_shannon_length) * 100; 
  fprintf('Efficiency (Shannon): %.2f%%\n', efficiency_shannon); 
  disp('---------------------------------'); 

%% Display Huffman's code average length and efficiency

% Calculate the average code length for Huffman 
  average_huffman_length = 0; 
  for i = 1:length(probabilities) 
    average_huffman_length = average_huffman_length + probabilities(i) * length(codeword{i}); 
  end
  disp('---------------------------------'); 
  fprintf('Code average length (Huffman): %.4f bits/symbol\n', average_huffman_length); 
% Calculate Huffman efficiency (should not exceed 100%) 
  efficiency_huffman = (entropy / average_huffman_length) * 100; 
  fprintf('Efficiency (Huffman): %.2f%%\n', efficiency_huffman); 
  disp('---------------------------------'); 

%% Compare between the input and output images 
  if isequal(original_image, huffman_img) 
    disp('No Lost Data (Huffman)'); 
  else 
    disp('Some data is lost (Huffman)'); 
  end

  if isequal(original_image, shannon_img) 
    disp('No Lost Data (Shannon)'); 
  else 
    disp('Some data is lost (Shannon)'); 
  end

  
