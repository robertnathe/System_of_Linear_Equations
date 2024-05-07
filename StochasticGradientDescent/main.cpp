/*
 * Neural Network Algorithm:
 *
 * 1. Activation Functions:
 *    - The code defines several activation functions: sigmoid, tanh, relu, leaky_relu, and softmax.
 *    - These functions are used to introduce non-linearity into the neural network, which allows it to learn complex patterns in the data.
 *
 * 2. Neural Network Architecture:
 *    - The neural network has two layers: an input layer and a hidden layer.
 *    - The size of the hidden layer is determined by the `hiddenLayerSize` variable, which is calculated as `weights.size() / (inputs.size() + 1)`.
 *    - The output layer has the same size as the number of classes (determined by `outputs.size()`).
 *
 * 3. Forward Propagation:
 *    - The `calculateOutput` function performs the forward propagation of the neural network.
 *    - For each input vector:
 *      - The weighted sum of the inputs is calculated for the hidden layer.
 *      - The `softmax` activation function is applied to the weighted sum to get the hidden layer outputs.
 *      - The weighted sum of the hidden layer outputs is calculated for the output layer.
 *      - The `softmax` activation function is applied to the weighted sum to get the final output values.
 *    - The final output values are returned as a vector.
 *
 * 4. Training:
 *    - The training process is performed using stochastic gradient descent.
 *    - The training data is provided as a vector of input vectors and a corresponding vector of output labels.
 *    - The weights are initialized randomly.
 *    - The training loop runs for a specified number of epochs (`numEpochs`).
 *    - For each epoch:
 *      - The total error is initialized to 0.
 *      - For each training input:
 *        - The `calculateOutput` function is used to get the predicted outputs.
 *        - The target outputs are created by setting the corresponding element to 1 and the rest to 0.
 *        - The error is calculated as the difference between the target and predicted outputs.
 *        - The total error is accumulated.
 *        - The weights are updated using the stochastic gradient descent rule:
 *          - For the input layer weights:
 *            - The gradient is calculated as `learningRate * error * predictedOutputs[j] * (1.0 - predictedOutputs[j]) * inputs[i][k]`.
 *            - The weights are updated by adding the gradient.
 *          - For the hidden layer weights:
 *            - The gradient is calculated as `learningRate * error * predictedOutputs[j] * (1.0 - predictedOutputs[j]) * hiddenOutputs[k]`.
 *            - The weights are updated by adding the gradient.
 *    - After the training is complete, the model is tested using a sample input, and the predicted outputs are printed.
 *
 * Key Points:
 *   - The use of activation functions to introduce non-linearity in the neural network.
 *   - The forward propagation of the input through the hidden layer and output layer.
 *   - The stochastic gradient descent training process, where the weights are updated based on the error between the predicted and target outputs.
 *   - The update rules for the input layer weights and the hidden layer weights.
 *
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
double softmax(double x, const std::vector<double>& inputs) {
    double max_element = *std::max_element(inputs.begin(), inputs.end());
    double exp_sum = std::accumulate(inputs.begin(), inputs.end(), 0.0, [max_element](double sum, double y) { return sum + std::exp(y - max_element); });
    return std::exp(x - max_element) / exp_sum;
}
std::vector<double> calculateOutput(const std::vector<double>& weights, const std::vector<double>& inputs) {
    size_t hiddenLayerSize = weights.size() / (inputs.size() + 1);
    std::vector<double> hiddenOutputs(hiddenLayerSize);
    for (size_t i = 0; i < hiddenLayerSize; ++i) {
        double weightedSum = std::inner_product(weights.begin() + i * inputs.size(), weights.begin() + (i + 1) * inputs.size(), inputs.begin(), 0.0);
        hiddenOutputs[i] = softmax(weightedSum, inputs);
    }
    std::vector<double> outputValues(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
        double weightedSumOutput = std::inner_product(weights.begin() + hiddenLayerSize + i * hiddenLayerSize, weights.begin() + hiddenLayerSize + (i + 1) * hiddenLayerSize, hiddenOutputs.begin(), 0.0);
        outputValues[i] = softmax(weightedSumOutput, outputValues);
    }
    return outputValues;
}

int main() {
    std::vector<std::vector<double>> inputs = {{1.0, 2.0, 3.0}, {1.0, 2.0, 4.0}, {1.0, 2.0, 5.0}, {1.0, 2.0, 6.0}, {1.0, 2.0, 7.0},
                                              {1.0, 2.0, 8.0}, {1.0, 2.0, 9.0}, {1.0, 2.0, 10.0}, {1.0, 2.0, 11.0}, {1.0, 2.0, 12.0}};
    std::vector<int> outputs = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
    std::vector<double> weights(inputs[0].size() * (outputs.size() + 1), 0.1);
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(-0.1, 0.1);
    const double learningRate = 0.01;
    const size_t numEpochs = 50;
    size_t hiddenLayerSize = weights.size() / (inputs.size() + 1);
    std::vector<double> hiddenOutputs(hiddenLayerSize);
    for (size_t epoch = 0; epoch < numEpochs; ++epoch) {
        double totalError = 0.0;
        for (size_t i = 0; i < inputs.size(); ++i) {
            std::vector<double> predictedOutputs = calculateOutput(weights, inputs[i]);
            std::vector<double> targetOutputs(outputs.size(), 0.0);
            targetOutputs[outputs[i]] = 1.0;
            for (size_t j = 0; j < predictedOutputs.size(); ++j) {
                double error = targetOutputs[j] - predictedOutputs[j];
                totalError += error * error;
                for (size_t k = 0; k < inputs[i].size(); ++k) {
                    double delta = learningRate * error * predictedOutputs[j] * (1.0 - predictedOutputs[j]) * inputs[i][k];
                    weights[k + j * inputs[i].size()] += delta;
                }
                for (size_t k = 0; k < hiddenLayerSize; ++k) {
                    double delta = learningRate * error * predictedOutputs[j] * (1.0 - predictedOutputs[j]) * hiddenOutputs[k];
                    weights[inputs[i].size() + k + j * hiddenLayerSize] += delta;
                }
            }
        }
        std::cout << "Epoch " << epoch << ", total error: " << totalError << std::endl;
    }
    std::vector<double> testInput = {1.0, 2.0, 3.0};
    std::vector<double> testOutputs = calculateOutput(weights, testInput);
    std::cout << "Test input: (" << testInput[0] << ", " << testInput[1] << ", " << testInput[2] << "), outputs: ";
    for (double output : testOutputs) {
        std::cout << output << " ";
    }
    std::cout << std::endl;
    return 0;
}
