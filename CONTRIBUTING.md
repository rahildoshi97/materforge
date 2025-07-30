# Contributing to MaterForge

Thank you for your interest in contributing to MaterForge! 
We welcome contributions from the community and appreciate your help in making this project better.

## Getting Started

### Prerequisites
- Python 3.10 or higher
- Git
- Basic knowledge of Python and scientific computing

### Development Setup

1. Fork the repository on GitLab
2. Clone your fork locally:
```bash
git clone https://i10git.cs.fau.de/your-username/materforge.git
cd materforge
```
3. Create a virtual environment and install dependencies:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e .[dev]
```
4. Install pre-commit hooks:
```bash
pre-commit install
```

## Making Changes

### Branch Naming
- Use descriptive branch names: `feature/add-new-property`, `fix/temperature-calculation`, `docs/update-readme`
- Use prefixes: `feature/`, `fix/`, `docs/`, `refactor/`, `test/`

### Development Workflow
1. Create a new branch from `master`:
```bash
git checkout master
git pull origin master
git checkout -b feature/your-feature-name
```
2. Make your changes and test your changes thoroughly
3. Commit your changes with clear messages
4. Push to your fork and create a merge request

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific categories
pytest tests/unit/
pytest tests/integration/
pytest tests/performance/

# Run with coverage
pytest --cov=src/materforge
```

### Code Standards
Python Style
- Follow PEP 8 guidelines
- Use type hints
- Maximum line length: 120 characters
- Use descriptive variable names

### Quality Tools
Run before submitting changes:
```bash
# Format code
black src/ tests/
isort src/ tests/

# Check linting
flake8 src/ tests/
mypy src/materforge/
```

## Submitting Changes

### Merge Request Requirements
1. Ensure all tests pass
2. Code follows style guidelines
3. Update documentation if necessary
4. Create a merge request with:
    - Clear title and description
    - Reference any related issues
    - List of changes made
    - Testing performed

### Merge Request Template
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement

## Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Manual testing performed
```

Thank you for contributing to MaterForge! 🚀
