using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.Statistics;
using NAudio.Wave;
using System.Numerics;

const string wavFilePath = "../../../../data/input.wav";
const double fLow = 400;       // нижняя граница частотного диапазона
const double fHigh = 3000;     // верхняя граница частотного диапазона
const double windowDuration = 0.05;  // секунд
const double hopDuration = 0.025;    // секунд
const double peakThreshold = 1.5;    // множитель для динамического порога

// 1. Загрузить WAV и подготовить сигнал
var (sampleRate, samples) = LoadWav(wavFilePath);

// 2. Параметры окон
int windowSizeSamples = (int)Math.Round(windowDuration * sampleRate);
int hopSizeSamples = (int)Math.Round(hopDuration * sampleRate);

if (windowSizeSamples <= 0 || hopSizeSamples <= 0)
{
    Console.WriteLine("Неверные параметры окна/шага.");
    return;
}

// 3. Подготовить список окон (фрагментов сигнала)
var windows = SliceSignal(samples, windowSizeSamples, hopSizeSamples);
int numWindows = windows.Count;
if (numWindows == 0)
{
    Console.WriteLine("Сигнал слишком короткий для выбранных параметров окна.");
    return;
}

// 4. Длина FFT — ближайшая степень двойки
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

// 6. Матрица 0/1: строки = окна, столбцы = частотные бинЫ
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

// 8. Найти столбцы, состоящие целиком из 1
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

// 9. Определить частоту сигнала
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

// Вывод результатов
Console.WriteLine($"Найденная частота сигнала: {signalFrequency:F2} Гц");
Console.WriteLine($"Количество окон: {numWindows}, количество бинов в диапазоне: {numBins}");
Console.WriteLine("Вектор столбцов (1 — устойчивый сигнал):");
Console.WriteLine(string.Join("", binaryColumnsVector.Select(v => v.ToString())));
Console.WriteLine("Матрица 0/1 (каждая строка — окно):");
var to = File.CreateText("out.txt");
for (int w = 0; w < numWindows; w++)
{
    var row = new char[numBins];
    for (int j = 0; j < numBins; j++)
    {
        row[j] = binaryMatrix[w, j] == 1 ? '1' : '0';
    }
    Console.WriteLine(new string(row));
    to.WriteLine(new string(row));
}

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
    // Выполнить FFT
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
