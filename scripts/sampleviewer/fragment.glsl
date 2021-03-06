#version 400

in vec3 fragmentColor;
in float lightIntensity;

out vec3 color;

void main() {
	vec3 ambientColor = 0.3 * fragmentColor;
	vec3 diffuseColor = lightIntensity * fragmentColor;
	color = ambientColor + diffuseColor;
}
