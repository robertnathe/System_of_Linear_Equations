/*
* This C++ program implements a simple feedforward neural network with one hidden layer.
* It demonstrates fundamental concepts of neural networks, including:
*   - Weight initialization using He initialization (suitable for ReLU activation)
*   - Forward propagation to calculate network outputs
*   - Backpropagation for training the network
*   - Activation functions:
*   - ReLU (Rectified Linear Unit) for hidden neurons
*   - Sigmoid for output neurons (followed by softmax for multi-class classification)
* This program is an Artificial Neural Network (ANN).
*   - Feedforward structure: data progresses from input to output layers without loops.
*   - No feedback connections: the output of a layer doesn't feed back to itself or previous layers.
*   - Statelessness: the network computes the output solely based on the current input, not past inputs.
*   - Layer structure and activation: fully connected layers with ReLU for hidden and sigmoid for output.
*   - Data handling: processes data in batches or individually without dependency on previous inputs.
* The program performs the following tasks:
*   1. Initializes the network with weights using He initialization.
*   2. Trains the network using backpropagation to adjust weights based on errors and gradients.
*   3. Tests the network on new data to evaluate its performance.
* Training might require tuning hyperparameters (learning rate, initialization, network architecture) for better results. 
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <numeric>
// Utility function to safely create a vector
template <typename T>
std::vector<T> safe_vector(size_t size) {
    if (size > std::numeric_limits<typename std::vector<T>::size_type>::max() / sizeof(T)) {
        throw std::length_error("Vector size exceeds safe limits");
    }
    return std::vector<T>(size);
}
// Activation functions and their derivatives
std::vector<double> softmax(const std::vector<double> &z) {
    std::vector<double> sm(z.size());
    double max_element = *std::max_element(z.begin(), z.end());
    double sum = 0.0;
    for (double score : z) sum += std::exp(score - max_element);
    for (size_t i = 0; i < z.size(); ++i) sm[i] = std::exp(z[i] - max_element) / sum;
    return sm;
}
double relu(double x) {
    return std::max(0.0, x);
}
double relu_derivative(double x) {
    return x > 0 ? 1.0 : 0.0;
}
double sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}
double sigmoid_derivative(double x) {
    return x * (1.0 - x);
}
// Initialize weights using He Normal Initialization
void initialize_weights(std::vector<double> &weights, size_t fan_in) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(0, std::sqrt(2.0 / fan_in));
    std::generate(weights.begin(), weights.end(), [&]() { return dis(gen); });
}
// Forward pass function to calculate network outputs
std::vector<double> calculate_output(const std::vector<double> &weights, const std::vector<double> &inputs, size_t hidden_layer_size, size_t output_size) {
    size_t input_size = inputs.size();
    std::vector<double> hidden_outputs(hidden_layer_size);
    std::vector<double> output_values(output_size);
    size_t w_index = 0;
    // Calculate hidden layer outputs using ReLU activation
    for (size_t i = 0; i < hidden_layer_size; ++i) {
        hidden_outputs[i] = 0.0;
        for (size_t j = 0; j < input_size; ++j) {
            hidden_outputs[i] += weights[w_index++] * inputs[j];
        }
        hidden_outputs[i] = relu(hidden_outputs[i]);
    }
    // Calculate output layer values using sigmoid activation
    size_t output_weight_start = input_size * hidden_layer_size;
    double bias = weights[weights.size() - 1];
    for (size_t i = 0; i < output_size; ++i) {
        output_values[i] = bias;
        for (size_t j = 0; j < hidden_layer_size; ++j) {
            output_values[i] += weights[output_weight_start++] * hidden_outputs[j];
        }
        output_values[i] = sigmoid(output_values[i]);
    }
    // Apply softmax to the output layer for multi-class classification
    return softmax(output_values);
}
// Prepare training data
void prepare_data(std::vector<std::vector<double>> &inputs, std::vector<int> &outputs) {
    inputs = {
        {1.0 / 12.0, 2.0 / 12.0, 3.0 / 12.0}, {1.0 / 12.0, 2.0 / 12.0, 4.0 / 12.0},
        {1.0 / 12.0, 2.0 / 12.0, 5.0 / 12.0}, {1.0 / 12.0, 2.0 / 12.0, 6.0 / 12.0},
        {1.0 / 12.0, 2.0 / 12.0, 7.0 / 12.0}, {1.0 / 12.0, 2.0 / 12.0, 8.0 / 12.0},
        {1.0 / 12.0, 2.0 / 12.0, 9.0 / 12.0}, {1.0 / 12.0, 2.0 / 12.0, 10.0 / 12.0},
        {1.0 / 12.0, 2.0 / 12.0, 11.0 / 12.0}, {1.0 / 12.0, 2.0 / 12.0, 1.0}
    };
    outputs = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
}
// Initialize the network by setting up the weights
void initialize_network(size_t input_size, size_t hidden_layer_size, size_t output_size, std::vector<double> &weights) {
    size_t total_weights = input_size * hidden_layer_size + hidden_layer_size * output_size + 1;
    weights = safe_vector<double>(total_weights);
    initialize_weights(weights, input_size);
}
// Training function using backpropagation
void train_network(const std::vector<std::vector<double>> &inputs, const std::vector<int> &outputs,
                   size_t hidden_layer_size, size_t output_size, std::vector<double> &weights,
                   size_t num_epochs, double learning_rate) {
    size_t input_size = inputs[0].size();
    for (size_t epoch = 0; epoch < num_epochs; ++epoch) {
        double total_error = 0.0;
        // Iterate over each training instance
        for (size_t i = 0; i < inputs.size(); ++i) {
            std::vector<double> hidden_outputs(hidden_layer_size);
            std::vector<double> predicted_outputs = calculate_output(weights, inputs[i], hidden_layer_size, output_size);
            std::vector<double> target_outputs(output_size, 0.0);
            target_outputs[outputs[i]] = 1.0;
            std::vector<double> output_deltas(output_size);
            for (size_t j = 0; j < output_size; ++j) {
                double error = target_outputs[j] - predicted_outputs[j];
                output_deltas[j] = error * sigmoid_derivative(predicted_outputs[j]);
                total_error += error * error;
            }
            // Calculate errors and deltas for each hidden neuron
            std::vector<double> hidden_deltas(hidden_layer_size);
            size_t hidden_weight_start = 0;
            size_t output_weight_start = input_size * hidden_layer_size;
            for (size_t h = 0; h < hidden_layer_size; ++h) {
                double hidden_error = 0.0;
                for (size_t j = 0; j < output_size; ++j) {
                    hidden_error += output_deltas[j] * weights[output_weight_start + j * hidden_layer_size + h];
                }
                hidden_deltas[h] = hidden_error * relu_derivative(hidden_outputs[h]);
            }
            // Update weights for the output layer
            for (size_t j = 0; j < output_size; ++j) {
                for (size_t h = 0; h < hidden_layer_size; ++h) {
                    weights[output_weight_start + j * hidden_layer_size + h] += learning_rate * output_deltas[j] * hidden_outputs[h];
                }
            }
            // Update weights for the hidden layer
            for (size_t h = 0; h < hidden_layer_size; ++h) {
                for (size_t k = 0; k < input_size; ++k) {
                    weights[hidden_weight_start + h * input_size + k] += learning_rate * hidden_deltas[h] * inputs[i][k];
                }
            }
            // Update bias weight
            weights.back() += learning_rate * std::accumulate(output_deltas.begin(), output_deltas.end(), 0.0);
        }
        // Report progress every 10 epochs
        if (epoch % 10 == 0) {
            std::cout << "Epoch " << epoch << ", Average Error: " << total_error / inputs.size() << std::endl;
        }
    }
}
// Testing function to evaluate the model
void test_network(const std::vector<double> &test_input, size_t hidden_layer_size, size_t output_size, const std::vector<double> &weights) {
    std::vector<double> test_outputs = calculate_output(weights, test_input, hidden_layer_size, output_size);
    std::cout << "Test Input: ";
    for (auto value : test_input) std::cout << value << " ";
    std::cout << "\nOutput Probabilities: ";
    for (auto value : test_outputs) std::cout << value << " ";
    std::cout << std::endl;
}
int main() {
    const size_t inputSize = 3;
    const size_t hiddenLayerSize = 20;
    const size_t outputSize = 3;
    const size_t numEpochs = 300;
    const double learningRate = 0.2;
    std::vector<std::vector<double>> inputs;
    std::vector<int> outputs;
    prepare_data(inputs, outputs);
    std::vector<double> weights;
    initialize_network(inputSize, hiddenLayerSize, outputSize, weights);
    train_network(inputs, outputs, hiddenLayerSize, outputSize, weights, numEpochs, learningRate);
    std::vector<double> testInput = {1.0 / 12.0, 2.0 / 12.0, 3.0 / 12.0};
    test_network(testInput, hiddenLayerSize, outputSize, weights);
    return 0;
}
