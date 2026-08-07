using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.Statistics;
using NAudio.Wave;
using System.Numerics;
using Accord.Math.Wavelets;
using MatrixConvolution;

const string wavFilePath = "../../../../data/input5.wav";
const double fLow = 400;     
const double fHigh = 3000;     
const double windowDuration = 0.05;  
const double hopDuration = 0.025;   
const double peakThreshold = 1.5;    


var (sampleRate, samples) = LoadWav(wavFilePath);
samples = WaveletDenoise(samples);   


int windowSizeSamples = (int)Math.Round(windowDuration * sampleRate);
int hopSizeSamples = (int)Math.Round(hopDuration * sampleRate);

if (windowSizeSamples <= 0 || hopSizeSamples <= 0)
{
    Console.WriteLine("Неверные параметры окна/шага.");
    return;
}

var windows = SliceSignal(samples, windowSizeSamples, hopSizeSamples);
int numWindows = windows.Count;
if (numWindows == 0)
{
    Console.WriteLine("Сигнал слишком короткий для выбранных параметров окна.");
    return;
}

int nFft = NextPowerOfTwo(windowSizeSamples);
var freqAxis = ComputeFrequencyAxis(nFft, sampleRate);

// 5. Индексы частот в диапазоне [fLow, fHigh]
var freqIndices = Enumerable.Range(0, freqAxis.Length)
    .Where(i => freqAxis[i] >= fLow && freqAxis[i] <= fHigh)
    .ToArray();

if (freqIndices.Length == 0)
{
    Console.WriteLine("Нет частотных бинов в указанном диапазоне.");
    return;
}

int numBins = freqIndices.Length;
int[,] binaryMatrix = new int[numWindows, numBins];

for (int w = 0; w < numWindows; w++)
{
    double[] windowed = ApplyHannWindow(windows[w]);
    var rowBits = AnalyzeSpectrumToBinaryRow(windowed, nFft, freqIndices, peakThreshold);
    for (int j = 0; j < numBins; j++)
    {
        binaryMatrix[w, j] = rowBits[j];
    }
}

int[] binaryColumnsVector = new int[numBins];
for (int j = 0; j < numBins; j++)
{
    bool allOnes = true;
    for (int w = 0; w < numWindows; w++)
    {
        if (binaryMatrix[w, j] != 1)
        {
            allOnes = false;
            break;
        }
    }
    binaryColumnsVector[j] = allOnes ? 1 : 0;
}

var onesIndices = Enumerable.Range(0, numBins)
    .Where(j => binaryColumnsVector[j] == 1)
    .ToArray();

double signalFrequency = -1;
if (onesIndices.Length > 0)
{
    double sum = 0;
    foreach (int j in onesIndices)
    {
        sum += freqAxis[freqIndices[j]];
    }
    signalFrequency = sum / onesIndices.Length;
}

Console.WriteLine($"Найденная частота сигнала: {signalFrequency:F2} Гц");
ProcessMatrix(binaryMatrix);
Console.WriteLine($"Количество окон: {numWindows}, количество бинов в диапазоне: {numBins}");
Console.WriteLine("Вектор столбцов (1 — устойчивый сигнал):");
Console.WriteLine(string.Join("", binaryColumnsVector.Select(v => v.ToString())));
Console.WriteLine("Матрица 0/1 (каждая строка — окно):");
static (int SampleRate, double[] Samples) LoadWav(string path)
{
    using var reader = new WaveFileReader(path);
    int sampleRate = reader.WaveFormat.SampleRate;
    int channels = reader.WaveFormat.Channels;
    var sampleProvider = reader.ToSampleProvider();

    var samples = new List<double>();
    float[] buffer = new float[reader.WaveFormat.BlockAlign * 1024];

    int read;
    while ((read = sampleProvider.Read(buffer, 0, buffer.Length)) > 0)
    {
        for (int i = 0; i < read; i += channels)
        {
            double sum = 0;
            int availableChannels = Math.Min(channels, read - i);
            for (int ch = 0; ch < availableChannels; ch++)
            {
                sum += buffer[i + ch];
            }
            samples.Add(sum / availableChannels);
        }
    }

    return (sampleRate, samples.ToArray());
}

static List<double[]> SliceSignal(double[] samples, int windowSize, int hopSize)
{
    var result = new List<double[]>();
    for (int start = 0; start + windowSize <= samples.Length; start += hopSize)
    {
        double[] window = new double[windowSize];
        Array.Copy(samples, start, window, 0, windowSize);
        result.Add(window);
    }
    return result;
}

static int NextPowerOfTwo(int value)
{
    int n = 1;
    while (n < value)
    {
        n <<= 1;
    }
    return n;
}

static double[] ComputeFrequencyAxis(int nFft, int sampleRate)
{
    int half = nFft / 2;
    var freqs = new double[half + 1];
    for (int i = 0; i <= half; i++)
    {
        freqs[i] = (double)i * sampleRate / nFft;
    }
    return freqs;
}

static double[] ApplyHannWindow(double[] data)
{
    int n = data.Length;
    var windowed = new double[n];
    for (int i = 0; i < n; i++)
    {
        double weight = 0.5 * (1 - Math.Cos(2 * Math.PI * i / (n - 1)));
        windowed[i] = data[i] * weight;
    }
    return windowed;
}

static int[] AnalyzeSpectrumToBinaryRow(double[] windowed, int nFft, int[] freqIndices, double thresholdFactor)
{
    
    var complex = new Complex[nFft];
    for (int i = 0; i < windowed.Length; i++)
    {
        complex[i] = new Complex(windowed[i], 0);
    }
    Fourier.Forward(complex, FourierOptions.Matlab);

    int half = nFft / 2;
    double[] magnitude = new double[half + 1];
    for (int i = 0; i <= half; i++)
    {
        magnitude[i] = complex[i].Magnitude;
    }

    int numBins = freqIndices.Length;
    var rowBits = new int[numBins];

    var localMagnitudes = freqIndices.Select(idx => magnitude[idx]).ToArray();
    double mean = localMagnitudes.Mean();
    double std = localMagnitudes.StandardDeviation();
    double dynamicThreshold = mean + thresholdFactor * std;

    for (int j = 0; j < numBins; j++)
    {
        int idx = freqIndices[j];
        double val = magnitude[idx];

        if (val < dynamicThreshold)
        {
            rowBits[j] = 0;
            continue;
        }

        double left = magnitude[Math.Max(idx - 1, 0)];
        double right = magnitude[Math.Min(idx + 1, magnitude.Length - 1)];
        rowBits[j] = (val >= left && val >= right) ? 1 : 0;
    }
    return rowBits;
}

static double[] WaveletDenoise(double[] signal)
{
    double[] filtered = new double[signal.Length];

    int window = 8; 

    for (int i = 0; i < signal.Length; i++)
    {
        double sum = 0;
        int count = 0;

        for (int j = -window; j <= window; j++)
        {
            int idx = i + j;
            if (idx >= 0 && idx < signal.Length)
            {
                sum += signal[idx];
                count++;
            }
        }

        filtered[i] = sum / count;
    }
    return filtered;
}



static void ProcessMatrix(int[,] matrix)
{
    string dataDir = "../../../../data/";

    Directory.CreateDirectory(
        Path.Combine(dataDir, "out5dwt"));

    int[] sampleSizes = { 40, 30, 20 };

    foreach (var sampleSize in sampleSizes)
    {
        ProcessSample(matrix, sampleSize, dataDir);
    }
}

static void ProcessSample(
    int[,] matrix,
    int sampleSize,
    string dataDir)
{
    var outputPath =
        Path.Combine(
            dataDir,
            "out5dwt",
            sampleSize.ToString());

    Directory.CreateDirectory(outputPath);

    Tools.SaveMatrixAsBmp(
        matrix,
        Path.Combine(outputPath,
        "source5dwt.bmp"));

    IConvolution alg =
        new ConvolutionFixedSize(
            Path.Combine(outputPath, "conv"));

    int[,] kernel =
    {
        {0,1,0},
        {0,1,0},
        {0,1,0}
    };

    alg.DoConvolution(
        matrix,
        kernel,
        out var convolutedMatrix);

    Tools.SaveMatrixAsBmp(
        convolutedMatrix,
        Path.Combine(outputPath,
        "convoluted5dwt.bmp"));

    for (int y = 0;
         y < matrix.GetLength(0) - sampleSize;
         y += sampleSize)
    {
        var sampleMatrix =
            new int[
                sampleSize,
                convolutedMatrix.GetLength(1)];

        int sampleMatrixSize =
            sampleSize *
            convolutedMatrix.GetLength(1);

        int sourceIndex =
            y *
            convolutedMatrix.GetLength(1);

        Array.Copy(
            convolutedMatrix,
            sourceIndex,
            sampleMatrix,
            0,
            sampleMatrixSize);

        var signals =
            Tools.DetectSignal(sampleMatrix);

        Tools.SaveMatrixAsBmp(
            sampleMatrix,
            signals,
            Path.Combine(
                outputPath,
                $"signals5dwt.{y}.one.bmp"));

        Tools.SaveMatrixAsBmp2(
            sampleMatrix,
            signals,
            Path.Combine(
                outputPath,
                $"signals5dwt.{y}.two.bmp"));

        File.AppendAllLines(
            Path.Combine(
                outputPath,
                "result5dwt.txt"),
            new[]
            {
                $"sample: {y} Signals: {signals.Count}"
            });
    }
}

   