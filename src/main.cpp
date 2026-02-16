#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>
#include <optional>

using namespace std;

#define obj_len 13

struct Object {
    int type;
    sf::Glsl::Vec4 color;
    int mat;
    sf::Glsl::Vec3 pos;
    sf::Glsl::Vec3 size;
    float extra;
};

// project path for Linux
std::string getExecutablePath() {
    std::string exePath = std::filesystem::read_symlink("/proc/self/exe").parent_path();
    // going up one level
    return std::filesystem::path(exePath).parent_path().string();
}

int main() {
    unsigned int w = 1920;
    unsigned int h = 1080;
    sf::Vector2u size(w, h);
    int mouseX = w / 2;
    int mouseY = h / 2;
    sf::Clock clock;
    sf::Shader shader_1;
    bool RT = false;
    float mouseSensitivity = 0.001f;
    float speed = 0.1;
    bool mouseHidden = true;
    int stillFrames = 1;
    sf::Vector2i localPosition;
    sf::Vector3f position = sf::Vector3f(0.0f, 0.0f, 1.0f);

    std::string exePath = getExecutablePath();

    using namespace sf::Glsl;
    Object objects[obj_len] = {
        { 1, Vec4(0.8, 0.8, 0.8, 0), 1, Vec3(5., -5., 2.5), Vec3(2.5, 2.5, 2.5), 0. },
        { 1, Vec4(0.85, 0.85, 0.85, 0.9), 3, Vec3(6., 14., 2.5), Vec3(2.5, 2.5, 2.5), 0. },
        { 2, Vec4(0.4, 0.4, 0.8, 0.0), 3, Vec3(9.5, 5.8, 1.5), Vec3(1.5, 1.5, 1.5), 0. },
        { 2, Vec4(0.95, 0.95, 0.95, 0.99), 3, Vec3(15.5, 7., 2.), Vec3(2., 17., 30.), 0. },
        { 2, Vec4(0.8, 0.4, 0.4, 0), 0, Vec3(-9., 0., 0.), Vec3(-9., 0., 5.), 1.5 },
        { 1, Vec4(0.8, 0.8, 0.8, -0.3), 2, Vec3(-3., 0., 2.5), Vec3(2.5, 2.5, 2.5), 0. },
        { 2, Vec4(0.7, 0.95, 0.7, 0.0), 3, Vec3(-8., 8., 1.), Vec3(1., 2., 1.5), 0. },
        { 1, Vec4(0.8, 0.8, 0.8, 0.05), 0, Vec3(1.65, 9.5, 1.), Vec3(1., 0.35, 0.35), 0. },
        { 2, Vec4(0.9, 0.9, 0.9, 0.9), 0, Vec3(9.5, 0.0, 1.5), Vec3(1.5, 1.5, 1.5), 0. },
        { 2, Vec4(0.95, 0.95, 0.95, 0.99), 3, Vec3(0., -10., 2.), Vec3(17., 2, 30.), 0. },
        { 2, Vec4(0.95, 0.95, 0.95, 0.99), 3, Vec3(0., 20., 2.), Vec3(17., 2, 30.), 0. },
    };

    sf::Texture skybox;
    if (!skybox.loadFromFile(std::string(exePath + "/assets/Skybox.jpg"))) {
        std::cerr << "Failed to load Skybox.jpg" << std::endl;
        return -1;
    }

    if (!shader_1.loadFromFile(std::string(exePath + "/shaders/Shader_1.vert"), std::string(exePath + "/shaders/Shader_1.frag")))
    {
        cout << "Error: shader not loaded";
    }

    sf::RenderWindow window(sf::VideoMode({ w, h }), "Ray Tracing", sf::Style::Titlebar | sf::Style::Close); //sf::State::Fullscreen
    window.setFramerateLimit(60);

    sf::Mouse::setPosition(sf::Vector2i(w / 2, h / 2), window);
    window.setMouseCursorVisible(false);

    sf::Vector2u windowSize = window.getSize();
    shader_1.setUniform("u_resolution", sf::Glsl::Vec2(windowSize.x, windowSize.y));
    shader_1.setUniform("u_skybox", skybox);

    sf::RenderTexture firstTexture(sf::Vector2(w, h));
    sf::Sprite firstTextureSprite = sf::Sprite(firstTexture.getTexture());
    sf::Sprite firstTextureSpriteFlipped = sf::Sprite(firstTexture.getTexture());
    firstTextureSpriteFlipped.setScale(sf::Vector2f(1, -1));
    firstTextureSpriteFlipped.setPosition(sf::Vector2f(0, h));

    sf::RenderTexture outputTexture(sf::Vector2(w, h));
    sf::Sprite outputTextureSprite = sf::Sprite(outputTexture.getTexture());
    sf::Sprite outputTextureSpriteFlipped = sf::Sprite(firstTexture.getTexture());
    outputTextureSpriteFlipped.setScale(sf::Vector2f(1, -1));
    outputTextureSpriteFlipped.setPosition(sf::Vector2f(0, h));

    sf::RectangleShape rect(sf::Vector2f(w, h));
    rect.setFillColor(sf::Color::Black);

    std::random_device rd;
    std::mt19937 e2(rd());

    while (window.isOpen())
    {
        bool wasdSS[6] = { false, false, false, false, false, false };
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
                window.close();
            else if (event->is<sf::Event::MouseMoved>())
            {
                if (mouseHidden)
                {
                    localPosition = sf::Mouse::getPosition(window);
                    int deltaX = localPosition.x - (w / 2);
                    int deltaY = localPosition.y - (h / 2);
                    mouseX += deltaX;
                    mouseY += deltaY;
                    if (deltaX != 0 || deltaY != 0) { stillFrames = 1; }
                }
            }
            else if (event->is<sf::Event::KeyPressed>())
            {
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Enter)) {
                    RT = RT ? false : true; stillFrames = 1;
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Escape)){
                    window.close();
                }
                // Movement
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::W)) { wasdSS[0] = true; stillFrames = 1; }
                else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S)) { wasdSS[2] = true; stillFrames = 1; }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A)) { wasdSS[1] = true; stillFrames = 1; }
                else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D)) { wasdSS[3] = true; stillFrames = 1; }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Space)) { wasdSS[4] = true; stillFrames = 1; }
                else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::LShift)) { wasdSS[5] = true; stillFrames = 1; }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Equal)) speed *= 10.;
                else if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Hyphen)) speed /= 10.;
            }
        }

        sf::Mouse::setPosition(sf::Vector2i(w / 2, h / 2), window);

        sf::Vector3f dir = sf::Vector3f(0.0f, 0.0f, 0.0f);
        sf::Vector3f dirTemp;
        if (wasdSS[0]) dir = sf::Vector3f(1.0f, 0.0f, 0.0f);
        else if (wasdSS[2]) dir = sf::Vector3f(-1.0f, 0.0f, 0.0f);
        if (wasdSS[1]) dir += sf::Vector3f(0.0f, -1.0f, 0.0f);
        else if (wasdSS[3]) dir += sf::Vector3f(0.0f, 1.0f, 0.0f);
        if (wasdSS[4]) dir += sf::Vector3f(0.0f, 0.0f, 1.0f);
        else if (wasdSS[5]) dir += sf::Vector3f(0.0f, 0.0f, -1.0f);

        dirTemp.z = dir.z * cos(-mouseY * mouseSensitivity) - dir.x * sin(mouseY * mouseSensitivity);
        dirTemp.x = dir.z * sin(mouseY * mouseSensitivity) + dir.x * cos(mouseY * mouseSensitivity);
        dirTemp.y = dir.y;
        dir.x = dirTemp.x * cos(-mouseX * mouseSensitivity) - dirTemp.y * sin(-mouseX * mouseSensitivity);
        dir.y = dirTemp.x * sin(-mouseX * mouseSensitivity) + dirTemp.y * cos(-mouseX * mouseSensitivity);
        dir.z = -dirTemp.z;
        position += dir * speed;

        float time = clock.getElapsedTime().asSeconds();

        float type[obj_len];
        sf::Glsl::Vec4 color[obj_len];
        float mat[obj_len];
        sf::Glsl::Vec3 pos[obj_len];
        sf::Glsl::Vec3 size[obj_len];
        float extra[obj_len];
        for (int i = 0; i < obj_len; i++) {
            type[i] = objects[i].type;
            color[i] = objects[i].color;
            mat[i] = objects[i].mat == 3 && !RT ? 0 : objects[i].mat;
            pos[i] = objects[i].pos;
            size[i] = objects[i].size;
            extra[i] = objects[i].extra;
        }

        shader_1.setUniform("u_sample_part", 1.0f / stillFrames);

        shader_1.setUniform("len", float(obj_len));
        shader_1.setUniformArray("type", type, obj_len);
        shader_1.setUniformArray("color", color, obj_len);
        shader_1.setUniformArray("mat", mat, obj_len);
        shader_1.setUniformArray("pos", pos, obj_len);
        shader_1.setUniformArray("size", size, obj_len);
        shader_1.setUniformArray("extra", extra, obj_len);

        std::uniform_real_distribution<> dist(0.0f, 1.0f);
        shader_1.setUniform("u_seed1", sf::Vector2f((float)dist(e2), (float)dist(e2)) * 999.0f);
        shader_1.setUniform("u_seed2", sf::Vector2f((float)dist(e2), (float)dist(e2)) * 999.0f);

        shader_1.setUniform("u_time", time);
        shader_1.setUniform("u_mouse", sf::Glsl::Vec2(mouseX * mouseSensitivity, mouseY * mouseSensitivity));
        shader_1.setUniform("u_pos", sf::Glsl::Vec3(position));

        shader_1.setUniform("RT", RT);

        if(stillFrames < 50){ // max sample parts limiter
            if (stillFrames % 2 == 1)
            {
                shader_1.setUniform("u_sample", firstTexture.getTexture());
                outputTexture.draw(firstTextureSpriteFlipped, &shader_1);
                window.draw(outputTextureSprite);
            }
            else
            {
                shader_1.setUniform("u_sample", outputTexture.getTexture());
                firstTexture.draw(outputTextureSpriteFlipped, &shader_1);
                window.draw(firstTextureSprite);
            }
            stillFrames++;
        }

        window.display();
    }
}