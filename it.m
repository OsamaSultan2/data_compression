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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and process the image
  img_path = selected_file;
  original_image = imread(img_path);

% Convert to grayscale if needed
  if size(original_image, 3) == 3
      original_image = rgb2gray(original_image);
  end
    
% Display original image
  figure('Name', 'Original Image', 'NumberTitle', 'off');
  imshow(original_image);
  title('Original Grayscale Image');

% Get image dimensions
  [rows, cols] = size(original_image);
  fprintf('Image dimensions: %d x %d\n', rows, cols);
  fprintf('Total pixels: %d\n', rows * cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

codeword_length = length(probabilities); 
codeword = cell(1, codeword_length); % Empty cell array for Huffman codewords 
x = zeros(codeword_length, codeword_length); % Matrix for Huffman algorithm 
% Temporary probability array to keep original values intact 
temp_p = probabilities; 

% Huffman algorithm 
for a = 1:codeword_length-1 
[~, idx] = sort(temp_p); % Sort probabilities 
% Combine the smallest two probabilities 
temp_p(idx(2)) = temp_p(idx(1)) + temp_p(idx(2));  
x(idx(1), a) = 11; % Assign value to the first smallest p 
x(idx(2), a) = 10; % Assign value to the second smallest p 
temp_p(idx(1)) = nan; % Mark the combined node 
end 
% Generate Huffman codewords based on the matrix 
for a = codeword_length-1:-1:1 
code_1 = find(x(:, a) == 11); 
code_0 = find(x(:, a) == 10); 
codeword{code_1} = [codeword{code_0}, '1']; 
codeword{code_0} = [codeword{code_0}, '0']; 
end 
% Display Huffman dictionary with ASCII code 
disp('---------------------------------'); 
disp('Huffman Coding: Symbol | Codeword'); 
for i = 1:length(symbols) 
fprintf('%d | %s\n', symbols(i), codeword{i}); 
end 
% Encode the data using Huffman codes 
encoded = strings(1, length(img_vector)); 
for i = 1:length(img_vector) 
for j = 1:length(symbols) 
if img_vector(i) == symbols(j) 
encoded(i) = codeword{j}; 
end 
end 
end 

% Transmitted Binary Code (Encoded Message - Huffman) 
str = strjoin(encoded); 
str = strrep(str, ' ', ''); 
disp('---------------------------------'); 
disp('Transmitted Binary Code (Encoded Message - Huffman):'); 
disp(str); 
% Decode the message using Huffman codes 
decoded = []; % Initialize decoded data array 
stream = []; % Initialize temporary stream 
for i = 1:length(encoded) 
stream = [stream, encoded(i)]; 
idx = find(strcmp(stream, codeword)); 
if ~isempty(idx) 
decoded = [decoded, symbols(idx)]; 
stream = []; 
end 
end 
disp('---------------------------------'); 
disp('Decoded Data (Decoded Message - Huffman):'); 
img = uint8(original_image);
huffman_img = reshape(uint8(decoded), size(img));

figure('Name', 'Huffman Decoded Image', 'NumberTitle', 'off');
imshow(huffman_img); 

%% Shannon Coding using alpha method 
% Sort the symbols and probabilities based on probabilities 
[sorted_probabilities, idx] = sort(probabilities, 'descend'); 
sorted_symbols = symbols(idx); 
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
decoded_shannon = []; 
stream_shannon = []; 
for i = 1:length(encoded_shannon) 
stream_shannon = [stream_shannon, encoded_shannon(i)]; 
idx_shannon = find(strcmp(stream_shannon, shannon_codeword)); 
if ~isempty(idx_shannon) 
decoded_shannon = [decoded_shannon, sorted_symbols(idx_shannon)]; 
stream_shannon = []; 
end 
end 
disp('---------------------------------'); 
disp('Decoded Data (Decoded Message - Shannon):'); 
%disp(decoded_shannon); 
% Calculate the average code length for Shannon 
average_shannon_length = sum(sorted_probabilities .* codeLengths);
shannon_img = reshape(uint8(decoded_shannon), size(img));

figure('Name', 'shannon Decoded Image', 'NumberTitle', 'off');
imshow(shannon_img);


%% Display Shannon's code average length and efficiency 
disp('---------------------------------'); 
fprintf('Code average length (Shannon): %.4f bits/symbol\n', average_shannon_length); 
% Calculate Shannon efficiency (should not exceed 100%) 
efficiency_shannon = (entropy / average_shannon_length) * 100; 
if efficiency_shannon > 100 
efficiency_shannon = 100; % Cap efficiency at 100% if necessary 
end 
fprintf('Efficiency (Shannon): %.2f%%\n', efficiency_shannon); 
disp('---------------------------------'); 
% Calculate the average code length for Huffman 
average_huffman_length = 0; 
for i = 1:length(probabilities) 
average_huffman_length = average_huffman_length + probabilities(i) * length(codeword{i}); 
end 
% Display Huffman's code average length and efficiency 
disp('---------------------------------'); 
fprintf('Code average length (Huffman): %.4f bits/symbol\n', average_huffman_length); 
% Calculate Huffman efficiency (should not exceed 100%) 
efficiency_huffman = (entropy / average_huffman_length) * 100; 
if efficiency_huffman > 100 
efficiency_huffman = 100; % Cap efficiency at 100% if necessary 
end 
fprintf('Efficiency (Huffman): %.2f%%\n', efficiency_huffman); 
disp('---------------------------------'); 
% Compare input and output files 
if strcmp(img_vector, decoded) 
disp('No Lost Data (Huffman)'); 
else 
disp('Some data is lost (Huffman)'); 
end 
if strcmp(img_vector, decoded_shannon) 
disp('No Lost Data (Shannon)'); 
else 
disp('Some data is lost (Shannon)'); 
end

  
