class Audio
{
    int band;
    String name;
    SoundFile input;
    float[] spectrum;
    FFT fft;
    Amplitude amp;

    Audio(int band, String name)
    {
        this.band = band;
        this.name = name;

        spectrum = new float[band];

        input = new SoundFile(fluid_simulation.this, name);
        fft = new FFT(fluid_simulation.this, band);
        amp = new Amplitude(fluid_simulation.this);

        fft.input(input);
        amp.input(input);
    }

    void play()
    {
        input.play();
    }

    int getSpectrumID(int index)
    {
        return index%band;
    }

    float getFrequency(int id)
    {
        return spectrum[id];
    }
    float getFrequency(int id, float mult)
    {
        return spectrum[id]*mult;
    }

    float getAmplitude()
    {
        return amp.analyze();
    }
    float getAmplitude(float mult)
    {
        return amp.analyze()*mult;
    }

    void updateSpectrum()
    {
        fft.analyze(spectrum);
    }

} 