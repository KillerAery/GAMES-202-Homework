#include "denoiser.h"

Denoiser::Denoiser() : m_useTemportal(false) {}

void Denoiser::Reprojection(const FrameInfo &frameInfo) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    Matrix4x4 preWorldToScreen =
        m_preFrameInfo.m_matrix[m_preFrameInfo.m_matrix.size() - 1];
    Matrix4x4 preWorldToCamera =
        m_preFrameInfo.m_matrix[m_preFrameInfo.m_matrix.size() - 2];

#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Reproject
            int id = frameInfo.m_id(x, y);

            if (id < 0) {
                m_valid(x, y) = false;
                continue;
            }

            Matrix4x4 ReprojectMatrix =
                preWorldToScreen * m_preFrameInfo.m_matrix[id] *
                Inverse(frameInfo.m_matrix[id]);

            Float3 preScreenPos =
                ReprojectMatrix(frameInfo.m_position(x, y), Float3::EType::Point);

            if (preScreenPos.x < 0 || preScreenPos.x >= width || preScreenPos.y < 0 ||
                preScreenPos.y >= height) {
                m_valid(x, y) = false;
                continue;
            }

            int preid = m_preFrameInfo.m_id(preScreenPos.x, preScreenPos.y); 
            if (preid != id) {
                m_valid(x, y) = false;
                continue;
            }

            m_valid(x, y) = true;
            m_misc(x, y) = m_accColor(preScreenPos.x, preScreenPos.y);
        }
    }
    std::swap(m_misc, m_accColor);
}

void Denoiser::TemporalAccumulation(const Buffer2D<Float3> &curFilteredColor) {
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    int kernelRadius = 3;
#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Temporal clamp
            Float3 aver = Float3(0, 0, 0);
            int sampleNum = 0;
            for (int dy = -kernelRadius; dy <= kernelRadius; ++dy) {
                int j = y + dy;
                if (j < 0 || j >= height) continue;
                for (int dx = -kernelRadius; dx <= kernelRadius; ++dx) {
                    int i = x + dx;
                    if (i < 0 || i >= width) continue;
                    aver += curFilteredColor(i, j);
                    sampleNum++;
                }
            }
            aver /= sampleNum;

            Float3 variance = Float3(0, 0, 0);
            for (int dy = -kernelRadius; dy <= kernelRadius; ++dy) {
                int j = y + dy;
                if (j < 0 || j >= height) continue;
                for (int dx = -kernelRadius; dx <= kernelRadius; ++dx) {
                    int i = x + dx;
                    if (i < 0 || i >= width) continue;
                    Float3 delta = curFilteredColor(i, j) - aver;
                    variance += delta * delta;
                }
            }
            variance /= sampleNum;

            Float3 preColor = m_accColor(x, y);
            preColor = Clamp(preColor, 
                aver - variance * m_colorBoxK,
                aver + variance * m_colorBoxK);
            // TODO: Exponential moving average
            float alpha = m_valid(x, y) ? m_alpha : 1.0f;
            m_misc(x, y) = Lerp(preColor, curFilteredColor(x, y), alpha);
        }
    }
    std::swap(m_misc, m_accColor);
}

Buffer2D<Float3> Denoiser::Filter(const FrameInfo &frameInfo) {
    int height = frameInfo.m_beauty.m_height;
    int width = frameInfo.m_beauty.m_width;
    Buffer2D<Float3> filteredImage = CreateBuffer2D<Float3>(width, height);
    int kernelRadius = 16;
#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // TODO: Joint bilateral filter
            Float3 color = Float3(0.0f,0.0f,0.0f);
            float cumulate = 0.0f;

            for (int dy = -kernelRadius; dy <= kernelRadius; ++dy) {
                int j = y + dy;
                if (j < 0 || j >= height) continue;
                for (int dx = -kernelRadius; dx <= kernelRadius; ++dx) {
                    int i = x + dx;
                    if (i < 0 || i >= width) continue;

                    float Dij2 = dx * dx + dy * dy;

                    Float3 deltaColor = frameInfo.m_beauty(x, y) - frameInfo.m_beauty(i, j);
                    float DCiCj2 = Dot(deltaColor, deltaColor);

                    float Dnormal2 = 0.0f;
                    Dnormal2 = SafeAcos(Dot(frameInfo.m_normal(x, y), frameInfo.m_normal(i,j)));
                    Dnormal2 = std::min(1.0f,std::max(0.0f, Dnormal2));
                    Dnormal2 *= Dnormal2;

                    float Dplane2 = 0.0f;
                    Float3 deltaPosition = frameInfo.m_position(i, j) - frameInfo.m_position(x, y);
                    if (Dot(deltaPosition, deltaPosition) > 0.0001f) {
                        Dplane2 = Dot(frameInfo.m_normal(x, y), Normalize(deltaPosition));
                        Dplane2 *= Dplane2;
                    }

                    float J = exp(-Dij2 / (2.0f * m_sigmaCoord * m_sigmaCoord) -
                                  DCiCj2 / (2.0f * m_sigmaColor * m_sigmaColor) -
                                  Dnormal2 / (2.0f * m_sigmaNormal * m_sigmaNormal) -
                                  Dplane2 / (2.0f * m_sigmaPlane * m_sigmaPlane));
                    cumulate += J;
                    color += frameInfo.m_beauty(i, j) * J;
                }
            }

            filteredImage(x, y) = color / cumulate;
        }
    }
    return filteredImage;
}

void Denoiser::Init(const FrameInfo &frameInfo, const Buffer2D<Float3> &filteredColor) {
    m_accColor.Copy(filteredColor);
    int height = m_accColor.m_height;
    int width = m_accColor.m_width;
    m_misc = CreateBuffer2D<Float3>(width, height);
    m_valid = CreateBuffer2D<bool>(width, height);
    // pink room 所需要的 sigmacolor 调整
    if (width == 1280) {
        m_sigmaColor = 10.0f;
    }
}

void Denoiser::Maintain(const FrameInfo &frameInfo) { m_preFrameInfo = frameInfo; }

Buffer2D<Float3> Denoiser::ProcessFrame(const FrameInfo &frameInfo) {
    // Filter current frame
    Buffer2D<Float3> filteredColor;
    filteredColor = Filter(frameInfo);

    // Reproject previous frame color to current
    if (m_useTemportal) {
        Reprojection(frameInfo);
        TemporalAccumulation(filteredColor);
    } else {
        Init(frameInfo, filteredColor);
    }

    // Maintain
    Maintain(frameInfo);
    if (!m_useTemportal) {
        //m_useTemportal = true;
    }
    return m_accColor;
}
